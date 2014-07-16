#if defined(A01_3A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Shortwave Interface to Edwards-Slingo radiation scheme.
!
! Purpose:
!   This subroutine interface the Edwards-Slingo radiation scheme
!   in the shortwave.
!
! Method:
!   Principally, arrays are transferred to the appropriate formats.
!   Separate subroutines are called for each physical process.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SWRAD(IERR                                          &
!                       Mixing Ratios
     &   , H2O, CO2, O3, O2_MIX_RATIO                                   &
     &   , CO2_DIM1, CO2_DIM2, CO2_3D, L_CO2_3D                         &
     &   , L_use_stochem_CH4, CH4_stochem                               &
!                       Pressure Fields
     &   , PSTAR                                                        &
     &   , P_LAYER_BOUNDARIES                                           &
     &   , P_LAYER_CENTRES                                              &
!                       Temperatures
     &   , TAC                                                          &
!                       Options for treating clouds
     &   , L_GLOBAL_CLOUD_TOP, GLOBAL_CLOUD_TOP, L_microphysics         &
     &   , L_INHOM_CLOUD, INHOM_CLOUD, DP_CORR_STRAT, DP_CORR_CONV      &
     &   , SIN_LATITUDE                                                 &
!                       Stratiform Cloud Fields
     &   , L_CLOUD_WATER_PARTITION, L_PC2                               &
     &   , LCA_AREA, LCA_BULK, LCCWC1, LCCWC2                           &
!                       Convective Cloud Fields
     &   , cca, cccwp, ccw, lcbase, ccb, cct                            &
!                       Surface Fields
     &   , LAND_ALBEDO, L_MOSES_II, l_cable, L_CTILE                    &
     &   , L_USE_SPEC_SEA                                               &
     &   , LAND_ALB, SICE_ALB, FLANDG                                   &
     &   , OPEN_SEA_ALBEDO, ICE_FRACTION, LAND, LAND0P5, LYING_SNOW     &
!                       Solar Fields
     &   , COSZIN, LIT, LIST, SCS                                       &
!                       Aerosol Fields
     &   , L_CLIMAT_AEROSOL, L_CLIM_AERO_HGT, L_HadGEM1_Clim_Aero       &
     &   , BL_DEPTH, N_LEVELS_BL                                        &
     &   , L_USE_CLEARRH, L_USE_DUST, DUST_DIM1, DUST_DIM2              &
     &   , DUST_1, DUST_2, DUST_3, DUST_4, DUST_5, DUST_6               &
     &   , L_USE_BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2, BIOGENIC       &
     &   , L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT                     &
     &   , SULP_DIM1, SULP_DIM2                                         &
     &   , ACCUM_SULPHATE, AITKEN_SULPHATE, DISS_SULPHATE               &
     &   , L_VOLCTS, VOLCMASS                                           &
     &   , SEA_SALT_FILM, SEA_SALT_JET, L_USE_SEASALT_INDIRECT          &
     &   , L_USE_SEASALT_DIRECT, SALT_DIM_A, SALT_DIM_B                 &
     &   , L_USE_SOOT_DIRECT, SOOT_DIM1, SOOT_DIM2                      &
     &   , FRESH_SOOT, AGED_SOOT                                        &
     &   , L_USE_BMASS_DIRECT, BMASS_DIM1, BMASS_DIM2                   &
     &   , FRESH_BMASS, AGED_BMASS, CLOUD_BMASS, L_USE_BMASS_INDIRECT   &
     &   , L_USE_OCFF_DIRECT, OCFF_DIM1, OCFF_DIM2                      &
     &   , FRESH_OCFF, AGED_OCFF, CLOUD_OCFF, L_USE_OCFF_INDIRECT       &
     &   , L_USE_ARCL, ARCL_DIM1, ARCL_DIM2, N_ARCL_SPECIES             &
     &   , N_ARCL_COMPNTS, I_ARCL_COMPNTS, ARCL                         &
     &   , AERO_MESO, L_MURK_RAD, Ntot_land, Ntot_sea                   &
!                       Level of tropopause
     &   , TRINDX                                                       &
!                       Spectrum
#include "swsarg3a.h"
!                       Algorithmic options
#include "swcarg3a.h"
     &   , PTS, L_MOD_K_FLUX, L_Wenyi                                   &
!                       Diagnostics
     &   , SW_diag, row_list, col_list                                  &
!                       Physical Dimensions
     &   , NLIT, N_POINTS, NLEVS, N_LAYER, NCLDS, NWET, NOZONE          &
     &   , row_length, rows, NPD_FIELD, NPD_PROFILE, NPD_LAYER          &
     &   , NPD_COLUMN, N_CCA_LEV                                        &
!                       Output
     &   , SURF_DOWN_SW                                                 &
     &   , FLUX_BELOW_690NM_SURF, L_FLUX_BELOW_690NM_SURF               &
     &   , NETSW, TOP_ABSORPTION, SWSEA, SWOUT, dirpar_out, l_direct_par&
           ! Variables needed to calculate layer masses
     &   , rho_r2, r_rho_levels, r_theta_levels                         &
     &   , q, qcl, qcf, qcf2, qrain, qgraup                             &
     &   , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
     &   )
!

      USE solinc_data, ONLY: sol_bearing, lg_orog_corr, orog_corr,      &
     &                       lg_f_orog, f_orog, L_orog
      USE swrdiag_mod, ONLY: StrSWDiag
!
!
      USE rad_switches_mod, ONLY:                                       &
          lrad_ccrad
!dhb599 20110812: get SC (moved from downstairs)
#include "swsc.h"

      IMPLICIT NONE
!
!
!
!     COMDECKS INCLUDED
#include "c_r_cp.h"
#include "c_g.h"
!!dhb599 20110812: (moved upstairs)
!#include "swsc.h"
#include "c_pi.h"
!     INTERNAL DIMENSIONS OF THE CODE
#include "dimfix3a.h"
!     SPECTRAL REGIONS
#include "spcrg3a.h"
!     ANGULAR INTEGRATION
#include "angint3a.h"
!     TREATMENT OF SCATTERING
#include "sctmth3a.h"
!     OPTIONS TO THE CODE ALTERABLE IN THE UM.
#include "swopt3a.h"
!     OPTIONS TO THE CODE FIXED IN THE UM.
#include "swfix3a.h"
!     SOLVERS
#include "solver3a.h"
!     ERROR FLAGS
#include "error3a.h"
!     UNIT NUMBERS FOR PRINTED OUTPUT
#include "stdio3a.h"
!     Cloud Overlap schemes
#include "clschm3a.h"
!
!
!     DUMMY ARGUMENTS
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     DIMENSIONS OF ARRAYS:
      Integer, Intent(IN) :: row_length
!                              Length of rows on each domain
      Integer, Intent(IN) :: rows
!                              Number of rows in the domain
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE                                                  &
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER                                                    &
!             ARRAY SIZES FOR LAYERS
     &   , NPD_COLUMN
!             NUMBER OF COLUMNS PER POINT
!
!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_POINTS                                                     &
!             Total number of points, including unlit ones
     &   , NWET                                                         &
!             NUMBER OF WET LEVELS
     &   , NOZONE                                                       &
!             NUMBER OF LEVELS WITH OZONE
     &   , NLEVS                                                        &
!             Number of layers in the main model
     &   , N_LAYER                                                      &
!             Number of layers seen in the radiation scheme
     &   , NCLDS                                                        &
!             NUMBER OF CLOUDY LEVELS
     &   , N_LEVELS_BL                                                  &
!             Number of layers occupied by boundary-layer aerosol
!             if L_CLIM_AERO_HGT is false.
     &   , N_CCA_LEV
!             NUMBER OF CONVECTIVE CLOUD LEVELS
!
!     SPECTRAL DATA:
#include "swspdc3a.h"
!
!
!
!     GASEOUS MIXING RATIOS
      REAL                                                              &
                !, INTENT(IN)
     &     H2O(NPD_FIELD, NWET)                                         &
!             MASS MIXING RATIO OF WATER
     &   , CO2                                                          &
!             MASS MIXING RATIO OF CO2
     &   , O3(NPD_FIELD, NOZONE)                                        &
!             MASS MIXING RATIOS OF OZONE
     &   , O2_MIX_RATIO
!             MASS MIXING RATIO OF OXYGEN
!
!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL                                                              &
                !, INTENT(IN)
     &     PSTAR(NPD_FIELD)                                             &
!             SURFACE PRESSURES
     &   , P_LAYER_BOUNDARIES(NPD_FIELD,0:NLEVS)                        &
!            PRESSURE AT BOUNDARIES OF LAYERS
     &   , P_LAYER_CENTRES(NPD_FIELD,0:NLEVS)                           &
!            PRESSURE AT CENTRES OF LAYERS
     &   , TAC(NPD_FIELD, NLEVS)
!             TEMPERATURES AT CENTRES OF LAYERS
!
!     INCIDENT SOLAR RADIATION:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NLIT                                                         &
!             NUMBER OF LIT POINTS
     &   , LIST(NPD_FIELD)
!             LIST OF LIT POINTS
      REAL                                                              &
                !, INTENT(IN)
     &     COSZIN(NPD_FIELD)                                            &
!             COSINES OF ZENITH ANGLE
     &   , SCS                                                          &
!             SCALING OF SOLAR INCIDENT FIELD
     &   , LIT(NPD_FIELD)
!             FRACTION OF TIME POINT IS LIT
!


      Integer, intent(in) :: lcbase(npd_field)
         Real, intent(in) ::    ccw(npd_field, nwet)

!     MICROPHYSICAL FLAG:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_MICROPHYSICS                                               &
!             FLAG FOR PARAMETRIZED MICROPHYSICS
     &,    l_mcr_qcf2                                                   &
                          ! Use second ice category
     &,    l_mcr_qrain                                                  &
                          ! Use prognostic rain
     &,    l_mcr_qgraup                                                 &
                          ! Use graupel
     &,    l_mixing_ratio ! Use mixing ratios in layer mass calculation
!
!     OPTIONS FOR TREATING CLOUDS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_GLOBAL_CLOUD_TOP                                           &
!             FLAG TO USE A GLOBAL VALUE FOR THE TOPS OF CLOUDS
!             TO ENSURE REPRODUCIBLE RESULTS
     &   , L_INHOM_CLOUD
!             FLAG TO USE SCALING FACTORS FOR INHOMOGENEOUS CLOUD
      INTEGER                                                           &
                !, INTENT(IN)
     &     GLOBAL_CLOUD_TOP
!             GLOBAL TOPMOST CLOUDY LAYER
      REAL                                                              &
                !, INTENT(IN)
     &     INHOM_CLOUD(NPD_CLOUD_COMPONENT)                             &
!             SCALING FACTORS FOR INHOMOGENEOUS CLOUD
     &   , DP_CORR_STRAT                                                &
!             Decorrelation pressure scale for large scale cloud
     &   , DP_CORR_CONV                                                 &
!             Decorrelation pressure scale for convective cloud
     &   , SIN_LATITUDE(NPD_FIELD)       
!             Full latitude field for decorrelation overlap calculation
!
!
!     PROPERTIES OF STRATIFORM CLOUDS:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLOUD_WATER_PARTITION                                      &
!             FLAG TO USE PROGNOSTIC CLOUD ICE CONTENTS
     &   , L_PC2
!             Flag to use PC2 cloud scheme
      REAL                                                              &
                !, INTENT(IN)
     &     LCCWC1(NPD_FIELD, NCLDS+1/(NCLDS+1))                         &
!             NOMINAL LIQUID WATER CONTENTS
     &   , LCCWC2(NPD_FIELD, NCLDS+1/(NCLDS+1))                         &
!             NOMINAL ICE WATER CONTENTS
     &   , LCA_AREA(NPD_FIELD, NCLDS+1/(NCLDS+1))                       &
!             AREA FRACTIONS OF LAYER CLOUDS OUTSIDE CONVECTIVE TOWERS
     &   , LCA_BULK(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             BULK FRACTIONS OF LAYER CLOUDS OUTSIDE CONVECTIVE TOWERS
!
!     PROPERTIES OF CONVECTIVE CLOUDS:

      INTEGER                                                           &
                !, INTENT(IN)
     &     CCB(NPD_FIELD)                                               &
!             BASE OF CONVECTIVE CLOUD
     &   , CCT(NPD_FIELD)
!             TOP OF CONVECTIVE CLOUD
      REAL                                                              &
                !, INTENT(IN)
     &     CCCWP(NPD_FIELD)                                             &
!             WATER PATH OF CONVECTIVE CLOUD
     &   , CCA(NPD_FIELD,N_CCA_LEV)
!             FRACTION OF CONVECTIVE CLOUD

!     AEROSOLS:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLIMAT_AEROSOL                                             &
!             FLAG FOR CLIMATOLOGICAL AEROSOL
     &   , L_CLIM_AERO_HGT                                              &
!             Flag to use the depth of the boundary layer to set
!             the climatological aerosol
     &   , L_HadGEM1_Clim_Aero                                          &
!             Flag to use HadGEM1 setting for climatological aerosols
     &   , L_MURK_RAD                                                   &
!             FLAG FOR MESOSCALE MODEL AEROSOL
     &   , L_Wenyi                                                      &
!             FLAG FOR WENYI P & T SCALING     
     &   , L_USE_CLEARRH
!             RELATIVE HUMIDITY USED FOR HYGROSCOPIC AEROSOLS
!             (GRID-BOX MEAN IF FALSE, CLEAR-SKY MEAN IF TRUE)
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_USE_SULPC_DIRECT                                           &
!             FLAG TO USE SULPHUR CYCLE FOR DIRECT EFFECT
     &   , L_USE_SULPC_INDIRECT                                         &
!             FLAG TO USE SULPHUR CYCLE FOR INDIRECT EFFECT
     &   , L_USE_DUST                                                   &
                      ! use direct rad effect of mineral dust
     &   , L_USE_BIOGENIC                                               &
                          ! use direct effect of biogenic aerosol
     &   , L_USE_SOOT_DIRECT                                            &
                             ! USE DIRECT RAD. EFFECT OF SOOT AEROSOL
     &   , L_USE_BMASS_DIRECT                                           &
                              ! USE DIRECT RAD. EFFECT OF BIOMASS SMOKE
     &   , L_USE_BMASS_INDIRECT                                         &
                                ! USE INDIRECT EFFECT OF BIOMASS SMOKE
     &   , L_USE_OCFF_DIRECT                                            &
                              ! USE DIRECT RAD. EFFECT OF OCFF
     &   , L_USE_OCFF_INDIRECT                                          &
                                ! USE INDIRECT EFFECT OF OCFF
     &   , L_USE_SEASALT_INDIRECT                                       &
!             FLAG TO USE SEA-SALT FOR INDIRECT EFFECT
     &   , L_USE_SEASALT_DIRECT
!             FLAG TO USE SEA-SALT FOR DIRECT EFFECT
      INTEGER                                                           &
                !, INTENT(IN)
     &     SULP_DIM1,SULP_DIM2                                          &
!             DIMENSIONS FOR _SULPHATE ARRAYS, (P_FIELD,P_LEVELS or 1,1)
     &   , DUST_DIM1, DUST_DIM2                                         &
!        dimensions for mineral dust arrays (p_field,p_levels or 1,1)
     &   , BIOGENIC_DIM1, BIOGENIC_DIM2                                 &
!        dimensions for biogenic array (p_field,p_levels or 1,1)
     &   , SOOT_DIM1, SOOT_DIM2                                         &
!          DIMENSIONS FOR SOOT ARRAYS (P_FIELD,P_LEVELS or 1,1)
     &   , BMASS_DIM1, BMASS_DIM2                                       &
!          DIMENSIONS FOR BIOMASS ARRAYS (P_FIELD,P_LEVELS or 1,1)
     &   , OCFF_DIM1, OCFF_DIM2                                         &
!          DIMENSIONS FOR OCFF ARRAYS (P_FIELD,P_LEVELS or 1,1)
     &   , SALT_DIM_A, SALT_DIM_B                                       &
!             DIMENSIONS FOR SALT ARRAYS ON INPUT (SALT_DIM_A=P_FIELD
!             AND SALT_DIM_B=P_LEVELS, OR ELSE 1,1)
     &   , SALT_DIM_IND_A, SALT_DIM_IND_B                               &
!             DIMENSIONS FOR SEA-SALT ARRAYS PASSED DOWN TO
!             R2_SET_CLOUD_FIELD IF INDIRECT EFFECT REQUIRED.
     &   , SALT_DIM_DIR_A, SALT_DIM_DIR_B
!             DIMENSIONS FOR SEA-SALT ARRAYS PASSED DOWN TO
!             R2_SET_AEROSOL_FIELD IF DIRECT EFFECT REQUIRED.
      REAL                                                              &
                !, INTENT(IN)
     &     ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)                         &
!             MASS MIXING RATIO OF ACCUMULATION MODE AEROSOL
     &   , AITKEN_SULPHATE(SULP_DIM1, SULP_DIM2)                        &
!             MASS MIXING RATIO OF AITKEN MODE AEROSOL
     &   , VOLCMASS(NPD_FIELD)                                          &
!             Mass of stratospheric volcanic aerosol at each point
     &   , DISS_SULPHATE(SULP_DIM1, SULP_DIM2)                          &
!             MIXING RATIO OF DISSOLVED SULPHATE
     &   , DUST_1(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div1 dust
     &   , DUST_2(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div2 dust
     &   , DUST_3(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div3 dust
     &   , DUST_4(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div4 dust
     &   , DUST_5(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div5 dust
     &   , DUST_6(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div6 dust
     &   , BIOGENIC(BIOGENIC_DIM1, BIOGENIC_DIM2)                       &
!             Mixing ratio of biogenic aerosol
     &   , SEA_SALT_FILM(SALT_DIM_A, SALT_DIM_B)                        &
!             Number concentration of film-mode sea-salt aerosol
     &   , SEA_SALT_JET(SALT_DIM_A, SALT_DIM_B)                         &
!             Number concentration of jet-mode sea-salt aerosol
     &,FRESH_SOOT(SOOT_DIM1,SOOT_DIM2),AGED_SOOT(SOOT_DIM1,SOOT_DIM2)   &
!             SOOT MIXING RATIOS
     &   , FRESH_BMASS(BMASS_DIM1, BMASS_DIM2)                          &
!             Mass mixing ratio of fresh biomass smoke
     &   , AGED_BMASS(BMASS_DIM1, BMASS_DIM2)                           &
!             Mass mixing ratio of aged biomass smoke
     &   , CLOUD_BMASS(BMASS_DIM1, BMASS_DIM2)                          &
!             Mass mixing ratio of in-cloud biomass smoke
     &   , FRESH_OCFF(OCFF_DIM1, OCFF_DIM2)                             &
!             Mass mixing ratio of fresh fossil-fuel organic carbon aer
     &   , AGED_OCFF(OCFF_DIM1, OCFF_DIM2)                              &
!             Mass mixing ratio of aged fossil-fuel organic carbon aer
     &   , CLOUD_OCFF(OCFF_DIM1, OCFF_DIM2)                             &
!             Mass mixing ratio of in-cloud fossil-fuel org carbon aer
     &   , BL_DEPTH(NPD_FIELD)                                          &
!             Depth of the boundary layer
     &   , AERO_MESO(NPD_FIELD, NLEVS)
!             MIXING RATIO OF 'URBAN' AEROSOL OF MESOSCALE MODEL
!
      LOGICAL   L_VOLCTS

!     AEROSOL CLIMATOLOGY FOR NWP
#include "arcl_dim.h"

      ! Number of requested species within the climatology
      INTEGER N_ARCL_SPECIES
      
      ! Corresponding number of requested components
      INTEGER N_ARCL_COMPNTS
      
      ! Model switch for each species
      LOGICAL L_USE_ARCL(NPD_ARCL_SPECIES)
      
      ! Index of each component
      INTEGER I_ARCL_COMPNTS(NPD_ARCL_COMPNTS)
      
      ! Array dimensions
      INTEGER                                                           &
     &        ARCL_DIM1                                                 &
     &   ,    ARCL_DIM2
     
      ! Mass-mixing ratios 
      REAL                                                              &
     &        ARCL(ARCL_DIM1, ARCL_DIM2, N_ARCL_COMPNTS)
!
!     CARBON CYCLE:
      LOGICAL   L_CO2_3D    !  controls use of 3D co2 field
      INTEGER                                                           &
                !, INTENT(IN)
     &     CO2_DIM1, CO2_DIM2
!             DIMENSIONS FOR CO2 ARRAY, (P_FIELD,P_LEVELS or 1,1)
      REAL                                                              &
                !, INTENT(IN)
     &     CO2_3D(CO2_DIM1, CO2_DIM2)
!             MASS MIXING RATIO OF CARBON DIOXIDE
!     STOCHEM:
      REAL                                                              &
                !, INTENT(IN)
     &     CH4_stochem(NPD_FIELD, NLEVS)
!             Mass mixing ratio of CH4 from STOCHEM
      Logical, Intent(IN) :: L_use_stochem_CH4
!
!     Properties of the surface:
      Logical, Intent(IN) :: Land(NPD_FIELD)
!                             Land mask
      Logical, Intent(IN) :: Land0p5(NPD_FIELD)
!                             Land mask (TRUE if land fraction > 0.5)
      Logical, Intent(IN) :: L_flux_below_690nm_surf
!                             Flag to calculate flux at wavelengths
!                             shorter than 690 nm at the surface:
!                             This may be required as a diagnostic
!                             or for MOSES
      Logical, Intent(IN) :: L_Moses_II
!                             Surface fluxes are required for MOSES II
      Logical, Intent(IN) :: l_cable   
!                             Surface fluxes are required for CABLE    
      Logical, Intent(IN) :: L_Ctile
!                             Switch for coastal tiling
      Logical, Intent(IN) :: L_USE_SPEC_SEA
!                             Switch for spectrally dep. sea albedos
      REAL                                                              &
                !, INTENT(IN)
     &     ICE_FRACTION(NPD_FIELD)                                      &
!             FRACTION OF SEA ICE IN SEA PORTION OF GRID BOX
     &   , LAND_ALBEDO(NPD_FIELD,4)                                     &
!             MOSES II LAND SURFACE ALBEDO FIELDS
!             (*,1) - DIRECT BEAM VISIBLE
!             (*,2) - DIFFUSE VISIBLE
!             (*,3) - DIRECT BEAM NEAR-IR
!             (*,4) - DIFFUSE NEAR-IR
     &   , LAND_ALB(NPD_FIELD)                                          &
!             SURFACE ALBEDO OF LAND
     &   , SICE_ALB(NPD_FIELD)                                          &
!             SURFACE ALBEDO OF SEA-ICE
     &   , FLANDG(NPD_FIELD)                                            &
!             Land fraction in grid box
     &   , OPEN_SEA_ALBEDO(NPD_FIELD, 2)                                &
!             SURFACE ALBEDO FIELD OF OPEN SEA
!             (DIRECT AND DIFFUSE COMPONENTS)
     &   , LYING_SNOW(NPD_FIELD)
!             MASS LOADING OF LYING SNOW
!
      Real                                                              &
                !, Intent(IN)
     &     Ntot_land                                                    &
                               ! Number of droplets over land / m-3
     &   , Ntot_sea            ! Number of droplets over sea / m-3
!
!                       Level of tropopause
      INTEGER                                                           &
     &     TRINDX(NPD_FIELD)
!             THE LAYER BOUNDARY OF THE TROPOPAUSE
!
!     INCREMENT OF TIME:
      REAL                                                              &
                !, INTENT(IN)
     &     PTS
!             TIME INCREMENT
!
!     Use modulus of fluxes to remove negative effective extinctions
      LOGICAL, INTENT(IN) :: L_MOD_K_FLUX

!     Information for the calculation of layer masses
      Real, intent(in)::                                                &
     &  rho_r2(npd_field,nlevs)                                         &
                                ! Air density*radius of earth**2 / kg m-1
     &, r_rho_levels(npd_field,nlevs)                                   &
                                      ! Height of rho levels / m
     &, r_theta_levels(npd_field,0:nlevs)                               &
                                           ! Height of theta levels / m
     &, q(npd_field,nwet)                                               &
                                ! Water vapour mixing ratio / kg kg-1
     &, qcl(npd_field,nwet)                                             &
                                ! Liquid water mixing ratio / kg kg-1
     &, qcf(npd_field,nwet)                                             &
                                ! Ice mixing ratio / kg kg-1
     &, qcf2(npd_field,nwet)                                            &
                                ! Second ice category mr / kg kg-1
     &, qrain(npd_field,nwet)                                           &
                                ! Rain mixing ratio / kg kg-1
     &, qgraup(npd_field,nwet)  ! Graupel mixing ratio / kg kg-1
!
!
!     CALCULATED FLUXES:
      REAL                                                              &
                !, INTENT(OUT)
     &     SWOUT(NPD_FIELD, NLEVS+2)                                    &
!             NET DOWNWARD FLUXES
     &   , SWSEA(NPD_FIELD)                                             &
!             SEA-SURFACE COMPONENTS OF FLUX
!             WEIGHTED BY (OPEN SEA)/(TOTAL SEA) FRACTION
     &   , NETSW(NPD_FIELD)                                             &
!             NET ABSORBED SHORTWAVE RADIATION
     &   , FLUX_BELOW_690NM_SURF(NPD_FIELD)                             &
!             NET SURFACE FLUX BELOW 690NM (AT POINTS WHERE THERE
!             IS SEA-ICE THIS IS WEIGHTED BY THE FRACTION OF OPEN SEA.)
!             NB: ONLY USED FOR NON MOSESII RUNS.
     &   , dirpar_out(npd_field)                                        &
!             Direct component of surface PAR flux
     &   , SURF_DOWN_SW(NPD_FIELD, 4)                                   &
!             SURFACE DOWNWARD SHORTWAVE RADIATION COMPONENTS
!             (*,1) - DIRECT BEAM VISIBLE
!             (*,2) - DIFFUSE VISIBLE
!             (*,3) - DIRECT BEAM NEAR-IR
!             (*,4) - DIFFUSE NEAR-IR
     &   , TOP_ABSORPTION(NPD_FIELD)
!             Radiative absorption above the top of the atmosphere
!             as seen in the main model
!
      LOGICAL l_direct_par
!             True if direct PAR flux required
!
!
!     Diagnostics:
!
      Type (StrSWDiag) :: SW_diag
!
      Integer, Intent(IN) :: row_list(npd_field)
!                              List of row indices of lit points
      Integer, Intent(IN) :: col_list(npd_field)
!                              List of column indices of lit points
!
!
!
!
!     LOCAL VARIABLES.
!
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      LOGICAL                                                           &
     &     L_CLEAR
!             CALCULATE CLEAR-SKY FIELDS
!     FLAGS FOR PROCESSES ACTUALLY ENABLED.
      LOGICAL                                                           &
     &     L_RAYLEIGH                                                   &
!             LOCAL FLAG FOR RAYLEIGH SCATTERING
     &   , L_GAS                                                        &
!             LOCAL FLAG FOR GASEOUS ABSORPTION
     &   , L_CONTINUUM                                                  &
!             LOCAL FLAG FOR CONTINUUM ABSORPTION
     &   , L_DROP                                                       &
!             LOCAL FLAG FOR SCATTERING BY DROPLETS
     &   , L_AEROSOL                                                    &
!             LOCAL FLAG FOR SCATTERING BY AEROSOLS
     &   , L_AEROSOL_CCN                                                &
!             LOCAL FLAG TO USE AEROSOLS TO DETERMINE CCN
     &   , L_ICE
!             LOCAL FLAG FOR SCATTERING BY ICE CRYSTALS
      INTEGER                                                           &
     &     I_SOLVER_CLEAR                                               &
!             SOLVER FOR CLEAR-SKY FLUXES
     &   , I_GAS_OVERLAP(NPD_BAND_SW)                                   &
!             OVERLAPS IN EACH BAND
     &   , I_SCATTER_METHOD_BAND(NPD_BAND_SW)
!             Treatment of scattering in each band
!
!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL                                                              &
     &     LATITUDE_G(NPD_PROFILE)                                      &
!             GATHERED LATITUDE ARRAY
     &   , D_MASS(NPD_PROFILE, NPD_LAYER)                               &
!             MASS THICKNESSES OF LAYERS
     &   , P(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             PRESSURE FIELD
     &   , T(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             TEMPERATURE FIELD
     &   , GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES_SW)     &
!             MASS FRACTIONS OF GASES
     &   , NULLMMR
!             NULL MASS MIXING RATIO
      PARAMETER(                                                        &
     &     NULLMMR=0.0E+00                                              &
     &   )
!
      REAL :: layer_heat_capacity(npd_profile, npd_layer)
!             Specific heat capacity of layer * d_mass

!     CLOUDY PROPERTIES:
      INTEGER                                                           &
     &     N_CONDENSED                                                  &
!             NUMBER OF CONDENSED PHASES
     &   , TYPE_CONDENSED(NPD_CLOUD_COMPONENT)                          &
!             TYPES OF CONDENSED COMPONENTS
     &   , I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)                       &
!             PARAMETRIZATION SCHEMES FOR COMPONENTS
     &   , N_CLOUD_TOP_GLOBAL                                           &
!             INVERTED GLOBAL TOPMOST CLOUDY LAYER
     &   , I_CLOUD_TMP
!             Cloud Overlap used by R2_CLOUD_LEVEL_DIAG
      REAL                                                              &
     &     CONDENSED_PARAM_LIST(NPD_CLOUD_PARAMETER_SW                  &
     &        , NPD_CLOUD_COMPONENT, NPD_BAND_SW)                       &
!             PARAMETERS FOR CONDENSED PHASES
     &   , CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER                 &
     &        , NPD_CLOUD_COMPONENT)                                    &
!             CHARACTERISTIC DIMENSIONS OF CONDENSED SPECIES
     &   , CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                &
     &        , NPD_CLOUD_COMPONENT)                                    &
!             MASS FRACTIONS OF CONDENSED SPECIES
     &   , W_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             CLOUD AMOUNTS
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             FRACTIONS OF DIFFERENT TYPES OF CLOUD
     &   , TOT_CLOUD_COVER(NPD_PROFILE)                                 &
!             Total cloud cover
     &   , CONDENSED_MIN_DIM(NPD_CLOUD_COMPONENT)                       &
!             MINIMUM DIMENSIONS OF CONDENSED COMPONENTS
     &   , CONDENSED_MAX_DIM(NPD_CLOUD_COMPONENT)
!             MAXIMUM DIMENSIONS OF CONDENSED COMPONENTS
!
!     PROPERTIES OF AEROSOLS:
      REAL                                                              &
     &     AEROSOL_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                  &
     &        , NPD_AEROSOL_SPECIES_SW)
!             MIXING RATIOS OF AEROSOLS
!
!     SOLAR FIELDS:
      REAL                                                              &
     &     SEC_0(NPD_PROFILE)                                           &
!             SECANTS OF ZENITH ANGLE
     &   , SOLAR_INCIDENT_NORM(NPD_PROFILE)
!             NORMALLY INCIDENT SOLAR IRRADIANCE
!
!     SURFACE PROPERTIES:
      LOGICAL                                                           &
     &     LAND_G(NPD_PROFILE)                                          &
!             GATHERED LAND MASK
     &   , LAND0P5_G(NPD_PROFILE)
!             GATHERED LAND MASK (TRUE if land fraction >0.5)
      INTEGER                                                           &
     &     I_SURFACE(NPD_PROFILE)
!             TYPES OF SURFACE AT EACH POINT
      REAL                                                              &
     &     ALBEDO_FIELD_DIFF_GREY(NPD_PROFILE)                          &
!             DIFFUSE ALBEDO FIELD
     &   , ALBEDO_FIELD_DIR_GREY(NPD_PROFILE)                           &
!             DIRECT ALBEDO FIELD
     &   , ALBEDO_FIELD_DIFF(NPD_PROFILE, NPD_BAND_SW)                  &
!             DIFFUSE ALBEDO FIELD
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE, NPD_BAND_SW)                   &
!             DIRECT ALBEDO FIELD
     &   , EMISSIVITY_FIELD(NPD_PROFILE, NPD_BAND_SW)                   &
!             EMISSIVITY FIELD
     &   , ALBEDO_SEA_DIFF_G(NPD_PROFILE, NPD_BAND_SW)                  &
!             GATHERED DIFFUSE ALBEDO FOR OPEN SEA
     &   , ALBEDO_SEA_DIR_G(NPD_PROFILE, NPD_BAND_SW)                   &
!             GATHERED DIRECT ALBEDO FOR OPEN SEA
     &   , FLANDG_G(NPD_PROFILE)                                        &
!             Gathered Land Fraction
     &   , ICE_FRACTION_G(NPD_PROFILE)
!             Gathered Fractional sea-ice
!
!     FLUXES:
      REAL                                                              &
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DIRECT FLUX
     &   , FLUX_DIFFUSE(NPD_PROFILE, 0:NPD_LAYER)                       &
!             DIFFUSE DOWNWARD FLUX     
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                 &
!             CLEAR-SKY DIRECT FLUX
     &   , FLUX_NET(NPD_PROFILE, 0: NPD_LAYER)                          &
!             NET/DOWNWARD FLUX
     &   , FLUX_NET_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                    &
!             CLEAR-SKY NET/DOWNWARD TOTAL FLUX
     &   , FLUX_UP(NPD_PROFILE, 0: NPD_LAYER)                           &
!             UPWARD FLUX
     &   , FLUX_UP_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR-SKY UPWARD FLUX
!
!     ARRAYS FOR USE WITH DIAGNOSTICS:
      REAL                                                              &
     &     WEIGHT_690NM(NPD_BAND_SW)
!             WEIGHTS FOR EACH BAND FOR REGION BELOW 690 NM
!
!     SURFACE FLUXES FOR COUPLING OR DIAGNOSTIC USE
      REAL                                                              &
     &     SEA_FLUX_G(NPD_PROFILE)                                      &
!             NET DOWNWARD FLUX INTO SEA
     &   , SURFACE_DOWN_FLUX_G(NPD_PROFILE)                             &
!             DOWNWARD FLUX AT SURFACE
     &   , SURF_DOWN_CLR_G(NPD_PROFILE)                                 &
!             CLEAR-SKY DOWNWARD FLUX AT SURFACE
     &   , SURF_UP_CLR_G(NPD_PROFILE)                                   &
!             CLEAR-SKY UPWARD FLUX AT SURFACE
     &   , FLUX_BELOW_690NM_SURF_G(NPD_PROFILE)                         &
!             GATHERED GRID BOX MEAN DOWNWARD SURFACE FLUX BELOW 690NM
     &   , FL_SEA_BELOW_690NM_SURF_G(NPD_PROFILE)
!             GATHERED OPEN SEA NET SURFACE FLUX BELOW 690NM
!
!     TEMPORARY FIELDS ASSOCIATED WITH 690NM FLUX OVER MEAN SOLID SURF.
      REAL                                                              &
     &     FRACSOLID                                                    &
                     ! land + sea-ice fraction in grid-box.
!
!
!        LOCAL FIELDS FOR CLOUD DIAGNOSTICS:
     &   , CLOUD_EXTINCTION_G(NPD_PROFILE, NPD_LAYER)                   &
!             MEAN EXTINCTION COEFFICIENT IN CLOUDS WEIGHTED BY THE
!             CLOUD AMOUNT AND THE CLEAR-SKY FLUX (GATHERED ARRAY)
     &   , CLOUD_WEIGHT_EXTINCTION_G(NPD_PROFILE, NPD_LAYER)            &
!             WEIGHTING FACTOR FOR EXTINCTION IN CLOUDS: THE PRODUCT
!             OF THE CLOUD AMOUNT AND THE CLEAR-SKY DIRECT FLUX
!             (GATHERED ARRAY)
     &   , LS_CLOUD_EXTINCTION_G(NPD_PROFILE, NPD_LAYER)                &
!             MEAN EXTINCTION COEFFICIENT IN CLOUDS WEIGHTED BY THE
!             CLOUD AMOUNT AND THE CLEAR-SKY FLUX (GATHERED ARRAY)
     &   , LS_CLOUD_WEIGHT_EXTINCTION_G(NPD_PROFILE, NPD_LAYER)         &
!             WEIGHTING FACTOR FOR EXTINCTION IN CLOUDS: THE PRODUCT
!             OF THE CLOUD AMOUNT AND THE CLEAR-SKY DIRECT FLUX
!             (GATHERED ARRAY)
     &   , CNV_CLOUD_EXTINCTION_G(NPD_PROFILE, NPD_LAYER)               &
!             MEAN EXTINCTION COEFFICIENT IN CLOUDS WEIGHTED BY THE
!             CLOUD AMOUNT AND THE CLEAR-SKY FLUX (GATHERED ARRAY)
     &   , CNV_CLOUD_WEIGHT_EXTINCTION_G(NPD_PROFILE, NPD_LAYER)        &
!             WEIGHTING FACTOR FOR EXTINCTION IN CLOUDS: THE PRODUCT
!             OF THE CLOUD AMOUNT AND THE CLEAR-SKY DIRECT FLUX
!             (GATHERED ARRAY)
     &   , SURF_VIS_DIR_G(NPD_PROFILE)                                  &
!             GATHERED DOWNWARD SURFACE DIRECT BEAM VISIBLE FLUX
     &   , SURF_VIS_DIF_G(NPD_PROFILE)                                  &
!             GATHERED DOWNWARD SURFACE DIFFUSE VISIBLE FLUX
     &   , SURF_NIR_DIR_G(NPD_PROFILE)                                  &
!             GATHERED DOWNWARD SURFACE DIRECT BEAM NEAR-INFRARED FLUX
     &   , SURF_NIR_DIF_G(NPD_PROFILE)
!             GATHERED DOWNWARD SURFACE DIFFUSE NEAR-INFRARED FLUX
!
!     TEMPORARY FIELD ASSOCIATED WITH THE OROGRAPHY CORRECTION.
      REAL                                                              &
     &     SWOUT_TEMP(NPD_FIELD)
!
!     FIELDS REQUIRED FOR CALL TO RADIATION CODE BUT NOT USED
      INTEGER                                                           &
     &     N_ORDER_GAUSS                                                &
     &   , I_GAS
!
!     AUXILIARY VARIABLES:
      REAL                                                              &
     &     CPBYG                                                        &
!             SPECIFIC HEAT BY GRAVITY
     &   , WEIGHT_BAND(NPD_BAND_SW)
!             WEIGHTING FACTORS FOR BANDS
      PARAMETER(CPBYG=CP/G)
!
!     VARIABLES REQUIRED FOR COMPATIBILITY WITH SUBROUTINES:
      INTEGER                                                           &
     &     N_FRAC_SOL_POINT                                             &
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)
      REAL                                                              &
     &     DUMMY1D(1)                                                   &
     &    ,DUMMY2D(1,1)                                                 &
     &    ,DUMMY3D(1,1,1)
!
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     R2_SET_GAS_MIX_RATIO, R2_SET_THERMODYNAMIC                   &
     &   , R2_SET_AEROSOL_FIELD, R2_SET_CLOUD_FIELD                     &
     &   , R2_SET_CLOUD_PARAMETRIZATION                                 &
     &   , R2_SET_SURFACE_FIELD_SW                                      &
     &   , R2_COMPARE_PROC
!
!
!
!      write(6,*)'XXX R2_SWRAD: SC = ', SC
!
!     INITIALIZE THE ERROR FLAG FOR THE RADIATION CODE.
      IERR=I_NORMAL

      IF (lrad_ccrad) THEN
        CONDENSED_DIM_CHAR(:,:,:) = 0.0E+00
      END IF 
!
!     INITIALIZATIONS FOR DIAGNOSTICS DEPENDING ON BANDS
!
      IF ( L_FLUX_BELOW_690NM_SURF .OR. L_MOSES_II .OR. l_cable ) THEN
! DEPENDS ON: r2_set_690nm_weight
         CALL R2_SET_690NM_WEIGHT(N_BAND_SW                             &
     &      , L_PRESENT_SW                                              &
     &      , N_BAND_EXCLUDE_SW                                         &
     &      , INDEX_EXCLUDE_SW                                          &
     &      , WAVE_LENGTH_SHORT_SW                                      &
     &      , WAVE_LENGTH_LONG_SW                                       &
     &      , WEIGHT_690NM                                              &
     &      , NPD_BAND_SW, NPD_EXCLUDE_SW, NPD_TYPE_SW                  &
     &      )
      ENDIF
!
!     COMPARE PROCESSES IN THE SPECTRAL FILE WITH THOSE ENABLED IN
!     THE CODE. ALSO WARN ABOUT COLLISION BETWEEN INTERACTIVE AEROSOLS
!     AND THE NWP AEROSOL CLIMATOLOGY.
! DEPENDS ON: r2_compare_proc
      CALL R2_COMPARE_PROC(IERR, L_PRESENT_SW                           &
     &   , L_RAYLEIGH_SW, L_GAS_SW, L_CONTINUUM_SW                      &
     &   , L_DROP_SW, L_AEROSOL_SW, L_AEROSOL_CCN_SW, L_ICE_SW          &
     &   , L_USE_DUST, L_USE_BIOGENIC                                   &
     &   , L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT                     &
     &   , L_USE_SEASALT_DIRECT                                         &
     &   , L_USE_SOOT_DIRECT                                            &
     &   , L_USE_BMASS_DIRECT                                           &
     &   , L_USE_OCFF_DIRECT                                            &
     &   , L_CLIMAT_AEROSOL                                             &
     &   , L_USE_ARCL, N_ARCL_SPECIES                                   &
     &   , L_MURK_RAD                                                   &
     &   , L_RAYLEIGH, L_GAS, L_CONTINUUM                               &
     &   , L_DROP, L_AEROSOL, L_AEROSOL_CCN, L_ICE                      &
     &   , NPD_TYPE_SW                                                  &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!
#if defined(SCMA)
       print 12,LAND_ALBEDO,LAND_ALB
12     format(1x,'usealb',4f6.3,2x,f6.3)
#endif
!
!     SET THE PROPERTIES OF THE SURFACE
! DEPENDS ON: r2_set_surface_field_sw
      CALL R2_SET_SURFACE_FIELD_SW(                                     &
     &     N_BAND_SW                                                    &
     &   , NLIT, LIST                                                   &
     &   , I_SURFACE, I_SPEC_SURFACE_SW                                 &
     &   , L_SURFACE_SW                                                 &
     &   , L_MICROPHYSICS, L_MOSES_II, l_cable, L_CTILE                 &
!     &   , L_MICROPHYSICS, L_MOSES_II,  L_CTILE                 &
     &   , L_USE_SPEC_SEA                                               &
     &   , LAND, LAND0P5, OPEN_SEA_ALBEDO                               &
     &   , LAND_ALB, SICE_ALB                                           &
     &   , FLANDG, ICE_FRACTION                                         &
     &   , LAND_ALBEDO, WEIGHT_690NM                                    &
     &   , EMISSIVITY_FIELD, ALBEDO_FIELD_DIR, ALBEDO_FIELD_DIFF        &
     &   , LAND_G, LAND0P5_G, FLANDG_G, ICE_FRACTION_G                  &
     &   , ALBEDO_SEA_DIFF_G, ALBEDO_SEA_DIR_G                          &
     &   , NPD_FIELD, NPD_PROFILE, NPD_BAND_SW, NPD_SURFACE_SW          &
     &   )
!
!     SET THE MIXING RATIOS OF GASES.
! DEPENDS ON: r2_set_gas_mix_ratio
      CALL R2_SET_GAS_MIX_RATIO(IERR                                    &
     &   , NLIT, NLEVS, N_LAYER, NWET, NOZONE                           &
     &   , LIST, L_EXTRA_TOP_SW                                         &
     &   , N_ABSORB_SW, TYPE_ABSORB_SW                                  &
     &   , .FALSE., .FALSE., .FALSE., .FALSE., L_O2_SW                  &
     &   , .FALSE., .FALSE., .FALSE., .FALSE.                           &
     &   , H2O, CO2, O3, NULLMMR, NULLMMR, NULLMMR, NULLMMR             &
     &   , O2_MIX_RATIO                                                 &
     &   , NULLMMR, NULLMMR, NULLMMR, NULLMMR                           &
     &   , GAS_MIX_RATIO                                                &
     &   , CO2_DIM1, CO2_DIM2, CO2_3D, L_CO2_3D                         &
     &   , CH4_stochem, L_use_stochem_CH4                               &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_SPECIES_SW            &
! Note that ozone is passed directly and all the other fields are 
! ignored here, so for the number of greenhouse gases we supply zero.
     &   ,0, dummy3d                                                    &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!     SET THE THERMODYNAMIC PROPERTIES OF THE ATMOSPHERE.
! DEPENDS ON: r2_set_thermodynamic
      CALL R2_SET_THERMODYNAMIC(NLIT, NLEVS, N_LAYER, nwet, LIST        &
     &   , L_EXTRA_TOP_SW, .FALSE.                                      &
     &   , PSTAR, DUMMY1D, DUMMY1D, DUMMY1D                             &
     &   , P_LAYER_BOUNDARIES                                           &
     &   , P_LAYER_CENTRES                                              &
     &   , DUMMY2D                                                      &
     &   , DUMMY2D, TAC                                                 &
     &   , rho_r2, r_rho_levels, r_theta_levels                         &
     &   , q, qcl, qcf, qcf2, qrain, qgraup                             &
     &   , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
     &   , P(:,1:), T(:,1:), DUMMY2D, DUMMY1D, DUMMY1D, DUMMY1D, D_MASS &
     &   , layer_heat_capacity                                          &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER                            &
     &   )
!
!
!     SET SEA-SALT ARRAY DIMENSIONS.
      IF (L_USE_SEASALT_DIRECT) THEN
         SALT_DIM_DIR_A=SALT_DIM_A
         SALT_DIM_DIR_B=SALT_DIM_B
      ELSE
         SALT_DIM_DIR_A=1
         SALT_DIM_DIR_B=1
      ENDIF
!
!
!     SET THE MIXING RATIOS OF AEROSOLS.
      IF (L_AEROSOL.OR.L_AEROSOL_CCN) THEN
! DEPENDS ON: r2_set_aerosol_field
         CALL R2_SET_AEROSOL_FIELD(IERR                                 &
     &      , NLIT, NLEVS, N_LAYER, N_AEROSOL_SW, TYPE_AEROSOL_SW       &
     &      , LIST, L_EXTRA_TOP_SW                                      &
     &      , L_CLIMAT_AEROSOL, L_CLIM_AERO_HGT, L_HadGEM1_Clim_Aero    &
     &      , BL_DEPTH, T, N_LEVELS_BL, L_MURK_RAD, AERO_MESO           &
     &      , L_USE_DUST, DUST_DIM1, DUST_DIM2                          &
     &      , DUST_1, DUST_2, DUST_3, DUST_4, DUST_5, DUST_6            &
     &      , L_USE_BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2              &
     &      , BIOGENIC                                                  &
     &      , L_USE_SULPC_DIRECT                                        &
     &      , SULP_DIM1, SULP_DIM2                                      &
     &      , ACCUM_SULPHATE, AITKEN_SULPHATE                           &
     &      , L_VOLCTS, VOLCMASS                                        &
     &      , L_USE_SEASALT_DIRECT, SALT_DIM_DIR_A, SALT_DIM_DIR_B      &
     &      , SEA_SALT_FILM, SEA_SALT_JET, P                            &
     &      , L_USE_SOOT_DIRECT, SOOT_DIM1, SOOT_DIM2                   &
     &      , FRESH_SOOT, AGED_SOOT                                     &
     &      , L_USE_BMASS_DIRECT, BMASS_DIM1, BMASS_DIM2                &
     &      , FRESH_BMASS, AGED_BMASS                                   &
     &      , L_USE_OCFF_DIRECT, OCFF_DIM1, OCFF_DIM2                   &
     &      , FRESH_OCFF, AGED_OCFF                                     &
     &      , N_ARCL_SPECIES, N_ARCL_COMPNTS, I_ARCL_COMPNTS            &
     &      , L_USE_ARCL, ARCL_DIM1, ARCL_DIM2, ARCL                    &
     &      , LAND0P5, LYING_SNOW, PSTAR                                &
     &      , P_LAYER_BOUNDARIES, TRINDX                                &
     &      , AEROSOL_MIX_RATIO                                         &
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES_SW &
     &      )
      ENDIF
!
!
!     ASSIGN THE PROPERTIES OF CLOUDS.
!
!
!     Check the consistency of cloud diagnostics.
      IF (SW_diag%re_conv_flag) THEN
         IF (.NOT.SW_diag%wgt_conv_flag) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** Error: Microphysical diagnostics for convective'    &
     &         , 'cloud must include the cloud weighting.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF ( (SW_diag%re_strat_flag).OR.(SW_diag%lwp_strat_flag) ) THEN
         IF (.NOT.SW_diag%wgt_strat_flag) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** Error: Microphysical diagnostics for stratiform'    &
     &         , 'cloud must include the cloud weighting.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (SW_diag%L_cloud_extinction) THEN
         IF (.NOT.SW_diag%L_cloud_weight_extinction) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** Error: The cloud extinction'                        &
     &           , 'may be diagnosed only in conjunction'               &
     &           , 'with the corresponding weights.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (SW_diag%L_ls_cloud_extinction) THEN
         IF (.NOT.SW_diag%L_ls_cloud_weight_extinction) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** Error: The layer cloud extinction'                  &
     &           , 'may be diagnosed only in conjunction'               &
     &           , 'with the corresponding weights.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (SW_diag%L_cnv_cloud_extinction) THEN
         IF (.NOT.SW_diag%L_cnv_cloud_weight_extinction) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** Error: The conv. cloud extinction'                  &
     &           , 'may be diagnosed only in conjunction'               &
     &           , 'with the corresponding weights.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
!
! DEPENDS ON: r2_set_cloud_parametrization
      CALL R2_SET_CLOUD_PARAMETRIZATION(IERR, N_BAND_SW                 &
     &   , I_ST_WATER_SW, I_CNV_WATER_SW, I_ST_ICE_SW, I_CNV_ICE_SW     &
     &   , L_DROP_TYPE_SW                                               &
     &   , I_DROP_PARAMETRIZATION_SW                                    &
     &   , DROP_PARAMETER_LIST_SW                                       &
     &   , DROP_PARM_MIN_DIM_SW, DROP_PARM_MAX_DIM_SW                   &
     &   , L_ICE_TYPE_SW                                                &
     &   , I_ICE_PARAMETRIZATION_SW                                     &
     &   , ICE_PARAMETER_LIST_SW                                        &
     &   , ICE_PARM_MIN_DIM_SW, ICE_PARM_MAX_DIM_SW                     &
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST                      &
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM                         &
     &   , NPD_BAND_SW                                                  &
     &   , NPD_DROP_TYPE_SW, NPD_ICE_TYPE_SW, NPD_CLOUD_PARAMETER_SW    &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!
!     SET SEA-SALT ARRAY DIMENSIONS.
      IF (L_USE_SEASALT_INDIRECT) THEN
         SALT_DIM_IND_A=SALT_DIM_A
         SALT_DIM_IND_B=SALT_DIM_B
      ELSE
         SALT_DIM_IND_A=1
         SALT_DIM_IND_B=1
      ENDIF
!
!
! DEPENDS ON: r2_set_cloud_field
      CALL R2_SET_CLOUD_FIELD(NLIT, NLEVS, N_LAYER, NCLDS               &
     &   , LIST                                                         &
     &   , P, T, D_MASS                                                 &
     &   , CCB, CCT, CCA, CCCWP, CCW, LCBASE                            &
     &   , LCCWC1, LCCWC2, LCA_AREA, LCA_BULK                           &
     &   , L_PC2, L_MICROPHYSICS, L_AEROSOL_CCN                         &
     &   , SEA_SALT_FILM, SEA_SALT_JET                                  &
     &   , L_USE_SEASALT_INDIRECT, SALT_DIM_IND_A, SALT_DIM_IND_B       &
     &   , L_USE_BIOGENIC, BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2       &
     &   , SULP_DIM1, SULP_DIM2, ACCUM_SULPHATE, DISS_SULPHATE          &
     &   , AITKEN_SULPHATE, L_USE_BMASS_INDIRECT                        &
     &   , BMASS_DIM1, BMASS_DIM2, AGED_BMASS, CLOUD_BMASS              &
     &   , L_USE_OCFF_INDIRECT, OCFF_DIM1, OCFF_DIM2                    &
     &   , AGED_OCFF, CLOUD_OCFF                                        &
     &   , LYING_SNOW                                                   &
     &   , L_CLOUD_WATER_PARTITION, LAND0P5_G, FLANDG_G                 &
     &   , I_CLOUD_SW, I_CLOUD_REPRESENTATION_SW, I_CONDENSED_PARAM     &
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM                         &
     &   , N_CONDENSED, TYPE_CONDENSED                                  &
     &   , W_CLOUD, FRAC_CLOUD, L_LOCAL_CNV_PARTITION_SW                &
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                      &
     &   , SW_diag%RE_CONV, SW_diag%RE_CONV_FLAG                        &
     &   , SW_diag%RE_STRAT, SW_diag%RE_STRAT_FLAG                      &
     &   , SW_diag%WGT_CONV, SW_diag%WGT_CONV_FLAG                      &
     &   , SW_diag%WGT_STRAT, SW_diag%WGT_STRAT_FLAG                    &
     &   , SW_diag%LWP_STRAT, SW_diag%LWP_STRAT_FLAG                    &
     &   , SW_diag%NTOT_DIAG, SW_diag%NTOT_DIAG_FLAG                    &
     &   , SW_diag%STRAT_LWC_DIAG, SW_diag%STRAT_LWC_DIAG_FLAG          &
     &   , SW_diag%SO4_CCN_DIAG, SW_diag%SO4_CCN_DIAG_FLAG              &
     &   , SW_diag%COND_SAMP_WGT, SW_diag%COND_SAMP_WGT_FLAG            &
     &   , SW_diag%Nc_diag, SW_diag%Nc_diag_flag                        &
     &   , SW_diag%Nc_weight, SW_diag%Nc_weight_flag                    &
     &   , col_list, row_list, row_length, rows                         &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES_SW    &
     &   , N_CCA_LEV, Ntot_land, Ntot_sea                               &
     &   )

      SELECT CASE (I_CLOUD_SW)
      CASE (IP_cloud_part_corr)
         I_CLOUD_TMP=IP_CLOUD_MIX_MAX
      CASE (IP_cloud_part_corr_cnv)
         I_CLOUD_TMP=IP_CLOUD_TRIPLE
      CASE DEFAULT
         I_CLOUD_TMP=I_CLOUD_SW
      END SELECT

      IF (SW_diag%weighted_re_flag.AND.                                 &
     &    SW_diag%sum_weight_re_flag) THEN
! DEPENDS ON: r2_cloud_level_diag
         CALL R2_CLOUD_LEVEL_DIAG(IERR, NLIT, N_LAYER, NCLDS            &
     &      , LIST                                                      &
     &      , I_CLOUD_TMP, I_CLOUD_REPRESENTATION_SW                    &
     &      , T, W_CLOUD, FRAC_CLOUD, .TRUE.                            &
     &      , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                   &
     &      , SW_diag%weighted_re_flag, SW_diag%weighted_re             &
     &      , SW_diag%sum_weight_re                                     &
     &      , col_list, row_list, row_length, rows                      &
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER                         &
     &      )
         IF (IERR /= I_NORMAL) RETURN
      ENDIF
!
      IF (SW_diag%wgtd_warm_re_flag.AND.                                &
     &    SW_diag%sum_wgt_warm_re_flag) THEN
! DEPENDS ON: r2_cloud_level_diag
         CALL R2_CLOUD_LEVEL_DIAG(IERR, NLIT, N_LAYER, NCLDS            &
     &      , LIST                                                      &
     &      , I_CLOUD_TMP, I_CLOUD_REPRESENTATION_SW                    &
     &      , T, W_CLOUD, FRAC_CLOUD, .FALSE.                           &
     &      , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                   &
     &      , SW_diag%wgtd_warm_re_flag, SW_diag%weighted_warm_re       &
     &      , SW_diag%sum_weight_warm_re                                &
     &      , col_list, row_list, row_length, rows                      &
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER                         &
     &      )
         IF (IERR /= I_NORMAL) RETURN
      ENDIF
!
!
      IF (L_orog) ALLOCATE(lg_orog_corr(NLIT))
!
!     SET THE INCIDENT SOLAR FLUX.
      DO L=1, NLIT
         SOLAR_INCIDENT_NORM(L)=SCS*SC*LIT(LIST(L))
         SEC_0(L)=1.0E+00/COSZIN(LIST(L))

! Gather the orography correction factor into lit points.
         IF (L_orog) THEN
            lg_orog_corr(L)=orog_corr(col_list(L),row_list(L))
         ENDIF
!
      ENDDO
!
!
!     CHECK THAT A VALID NUMBER HAS BEEN SUPPLIED FOR THE SOLVER.
      IF ( (I_SOLVER_SW /= IP_SOLVER_PENTADIAGONAL).AND.                &
     &     (I_SOLVER_SW /= IP_SOLVER_MIX_11).AND.                       &
     &     (I_SOLVER_SW /= IP_SOLVER_MIX_DIRECT).AND.                   &
     &     (I_SOLVER_SW /= IP_SOLVER_MIX_DIRECT_HOGAN).AND.             &
     &     (I_SOLVER_SW /= IP_SOLVER_HOMOGEN_DIRECT).AND.               &
     &     (I_SOLVER_SW /= IP_SOLVER_TRIPLE).AND.                       &
     &     (I_SOLVER_SW /= IP_SOLVER_TRIPLE_HOGAN)                      &
     &   ) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: AN INVALID SOLVER HAS BEEN SELECTED '           &
     &      , 'IN THE SHORTWAVE REGION.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!

!     SET CLEAR-SKY CALCULATIONS.
      L_CLEAR=SW_diag%L_solar_out_clear.OR.                             &
     &        SW_diag%L_surf_down_clr.OR.                               &
     &        SW_diag%L_clear_hr.OR.                                    &
     &       (SW_diag%L_cloud_extinction.AND.                           &
     &              SW_diag%L_cloud_weight_extinction)
!
      IF (L_CLEAR) THEN
!
!        SELECT A CLEAR-SKY SOLVER TO MATCH THE MAIN SOLVER.
         IF (I_SOLVER_SW == IP_SOLVER_PENTADIAGONAL) THEN
            I_SOLVER_CLEAR=IP_SOLVER_PENTADIAGONAL
         ELSE IF (I_SOLVER_SW == IP_SOLVER_MIX_11) THEN
            I_SOLVER_CLEAR=IP_SOLVER_PENTADIAGONAL
         ELSE IF (I_SOLVER_SW == IP_SOLVER_MIX_DIRECT) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ELSE IF (I_SOLVER_SW == IP_SOLVER_MIX_DIRECT_HOGAN) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ELSE IF (I_SOLVER_SW == IP_SOLVER_HOMOGEN_DIRECT) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ELSE IF (I_SOLVER_SW == IP_SOLVER_TRIPLE) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ELSE IF (I_SOLVER_SW == IP_SOLVER_TRIPLE_HOGAN) THEN
            I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
         ENDIF
!
      ENDIF
!
!
!     SET PROPERTIES FOR INDIVIDUAL BANDS.
      DO I=1, N_BAND_SW
         WEIGHT_BAND(I)=1.0E+00
         I_GAS_OVERLAP(I)=I_GAS_OVERLAP_SW
      ENDDO
!
!     Set the treatment of scattering.
      DO I=1, N_BAND_SW
         I_SCATTER_METHOD_BAND(I)=I_SCATTER_METHOD_SW
      ENDDO
!
!
!
!     INVERT THE TOPMOST CLOUDY LAYER IF USING A GLOBAL VALUE.
      IF (L_GLOBAL_CLOUD_TOP) THEN
         N_CLOUD_TOP_GLOBAL=N_LAYER+1-GLOBAL_CLOUD_TOP
      ENDIF
!
!     GATHER THE LATITUDE FIELD.
      DO L=1, NLIT
        LATITUDE_G(L)=ASIN(SIN_LATITUDE(LIST(L)))*180/PI
      ENDDO
!
!
! DEPENDS ON: flux_calc
      CALL FLUX_CALC(IERR                                               &
!                       Logical Flags for Processes
     &   , L_RAYLEIGH, L_AEROSOL, L_GAS, L_CONTINUUM                    &
     &   , L_CLOUD_SW, L_DROP, L_ICE, L_PC2                             &
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION_SW, I_2STREAM_SW, L_2_STREAM_CORRECT_SW&
     &   , L_RESCALE_SW, N_ORDER_GAUSS                                  &
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_BAND                                        &
!                       Options for treating clouds
     &   , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL                       &
     &   , L_INHOM_CLOUD, INHOM_CLOUD                                   &
!    &   , L_INHOM_CLOUD, INHOM_CLOUD                                   &
!                       Options for Solver
     &   , I_SOLVER_SW                                                  &
!                       General Spectral Properties
     &   , N_BAND_SW, 1, N_BAND_SW, WEIGHT_BAND                         &
!                       General Atmospheric Properties
     &   , NLIT, N_LAYER                                                &
     &   , L_LAYER_SW, L_CLOUD_LAYER_SW                                 &
     &   , P, T, DUMMY1D, DUMMY1D, DUMMY1D, DUMMY2D, D_MASS             &
!                       Spectral Region
     &   , ISOLIR_SW                                                    &
!                       Solar Fields
     &   , SEC_0, SOLAR_INCIDENT_NORM, SOLAR_FLUX_BAND_SW               &
     &   , RAYLEIGH_COEFFICIENT_SW                                      &
!                       Infra-red Fields
     &   , N_DEG_FIT_SW                                                 &
     &   , THERMAL_COEFFICIENT_SW                                       &
     &   , T_REF_PLANCK_SW, .FALSE.                                     &
!                       Gaseous Absorption
     &   , N_ABSORB_SW, I_GAS_OVERLAP, I_GAS                            &
     &   , GAS_MIX_RATIO                                                &
     &   , N_BAND_ABSORB_SW, INDEX_ABSORB_SW                            &
     &   , I_BAND_ESFT_SW                                               &
     &   , W_ESFT_SW, K_ESFT_SW                                         &
     &   , I_SCALE_ESFT_SW, I_SCALE_FNC_SW                              &
     &   , L_WENYI, SCALE_VECTOR_SW                                     &
     &   , P_REFERENCE_SW, T_REFERENCE_SW, L_MOD_K_FLUX                 &
!                       Doppler Broadening
     &   , L_DOPPLER_PRESENT_SW                                         &
     &   , DOPPLER_CORRECTION_SW                                        &
!                       Surface Fields
     &   , L_SURFACE_SW, I_SURFACE                                      &
     &   , I_SPEC_SURFACE_SW                                            &
     &   , SURFACE_ALBEDO_SW                                            &
     &   , ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR                          &
     &   , N_DIR_ALBEDO_FIT_SW                                          &
     &   , DIRECT_ALBEDO_PARM_SW                                        &
     &   , EMISSIVITY_GROUND_SW                                         &
     &   , EMISSIVITY_FIELD                                             &
!                       Continuum Absorption
     &   , N_BAND_CONTINUUM_SW                                          &
     &   , INDEX_CONTINUUM_SW, INDEX_WATER_SW                           &
     &   , K_CONTINUUM_SW, I_SCALE_FNC_CONT_SW                          &
     &   , SCALE_CONTINUUM_SW                                           &
     &   , P_REF_CONTINUUM_SW                                           &
     &   , T_REF_CONTINUUM_SW                                           &
!                       Properties of Aerosols
     &   , N_AEROSOL_SW                                                 &
     &   , AEROSOL_MIX_RATIO                                            &
     &   , AEROSOL_ABSORPTION_SW                                        &
     &   , AEROSOL_SCATTERING_SW                                        &
     &   , AEROSOL_ASYMMETRY_SW                                         &
     &   , I_AEROSOL_PARAMETRIZATION_SW                                 &
     &   , NHUMIDITY_SW                                                 &
     &   , HUMIDITIES_SW                                                &
     &   , TYPE_AEROSOL_SW                                              &
     &   , L_USE_CLEARRH                                                &
!                       Aerosol optical depth
!                       (computed in the LW: dummy values for SW)
     &   , N_AOD_WAVEL_SW                                               &
     &   , .FALSE., DUMMY2D, .FALSE., DUMMY2D, .FALSE., DUMMY2D         &
     &   , .FALSE., DUMMY2D, .FALSE., DUMMY2D, .FALSE., DUMMY2D         &
     &   , .FALSE., DUMMY2D, .FALSE., DUMMY2D                           &
     &   , AOD_ABSORPTION_SW, AOD_SCATTERING_SW, I_AOD_TYPE_SW          &
!                       Properties of Clouds
     &   , N_CONDENSED, TYPE_CONDENSED                                  &
     &   , I_CLOUD_SW, I_CLOUD_REPRESENTATION_SW, W_CLOUD, FRAC_CLOUD   &
     &   , TOT_CLOUD_COVER                                              &
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                      &
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST                      &
     &   , DP_CORR_STRAT, DP_CORR_CONV                                  &
     &   , LATITUDE_G                                                   &
!                       Fluxes Calculated
     &   , FLUX_DIRECT, FLUX_DIFFUSE, FLUX_NET, FLUX_UP                 &
     &   , DUMMY2D, DUMMY2D, DUMMY2D                                    &
     &   , SW_diag%L_flux_diffuse, .FALSE., .FALSE., .FALSE.            &
!                       Options for Clear-sky Fluxes
     &   , L_CLEAR, I_SOLVER_CLEAR                                      &
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_NET_CLEAR, FLUX_UP_CLEAR             &
!                       Arrays specific to the UM
!                       Arrays for Coupling
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION_G           &
     &   , ALBEDO_SEA_DIFF_G, ALBEDO_SEA_DIR_G, FLANDG_G                &
     &   , SEA_FLUX_G                                                   &
!                       Arrays for diagnostics specific to the UM
     &   , L_FLUX_BELOW_690NM_SURF, WEIGHT_690NM, DUMMY1D               &
     &   , FLUX_BELOW_690NM_SURF_G, FL_SEA_BELOW_690NM_SURF_G           &
     &   , L_MOSES_II,  L_CTILE                                         &
!     &   , L_MOSES_II, l_cable, L_CTILE                                 &
     &   , SURF_VIS_DIR_G,SURF_VIS_DIF_G,SURF_NIR_DIR_G,SURF_NIR_DIF_G  &
     &   , SW_diag%L_surface_down_flux, SURFACE_DOWN_FLUX_G             &
     &   , SW_diag%L_surf_down_clr, SURF_DOWN_CLR_G                     &
     &   , SW_diag%L_surf_up_clr, SURF_UP_CLR_G                         &
     &   , SW_diag%L_cloud_extinction, CLOUD_EXTINCTION_G               &
     &   , CLOUD_WEIGHT_EXTINCTION_G                                    &
     &   , SW_diag%L_ls_cloud_extinction, LS_CLOUD_EXTINCTION_G         &
     &   , LS_CLOUD_WEIGHT_EXTINCTION_G                                 &
     &   , SW_diag%L_cnv_cloud_extinction, CNV_CLOUD_EXTINCTION_G       &
     &   , CNV_CLOUD_WEIGHT_EXTINCTION_G                                &
     &   , .FALSE., DUMMY2D, DUMMY2D, .FALSE., DUMMY2D, DUMMY2D         &
     &   , .FALSE., DUMMY2D, DUMMY2D                                    &
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN, NPD_FIELD                &
     &   , NPD_BAND_SW                                                  &
     &   , NPD_SPECIES_SW                                               &
     &   , NPD_ESFT_TERM_SW, NPD_SCALE_FNC_SW                           &
     &   , NPD_SCALE_VARIABLE_SW                                        &
     &   , NPD_CONTINUUM_SW                                             &
     &   , NPD_AEROSOL_SPECIES_SW                                       &
     &   , NPD_HUMIDITIES_SW                                            &
     &   , NPD_CLOUD_PARAMETER_SW                                       &
     &   , NPD_THERMAL_COEFF_SW                                         &
     &   , NPD_SURFACE_SW, NPD_ALBEDO_PARM_SW                           &
     &   , NPD_AOD_WAVEL_SW                                             &
     &   ,l_cable                                                       &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!
!     PREPARE THE OUTPUT ARRAYS:
!
!     The former zeroing of the output arrays was taken up into
!     RAD_CTL at 5.3.
!
!     SCATTER THE NET DOWNWARD FLUX AT EACH LEVEL INTO SWOUT.
      DO I=1, NLEVS+1
         DO L=1, NLIT
            SWOUT(LIST(L), I)=FLUX_NET(L, N_LAYER+1-I)
         ENDDO
      ENDDO
!
!
!     NET SHORTWAVE RADIATION ABSORBED BY THE PLANET
!     (I. E. EARTH AND ATMOSPHERE TOGETHER):
!
      DO L=1, NLIT
         NETSW(LIST(L))=SWOUT(LIST(L), NLEVS+1)
      ENDDO
!
      IF (L_EXTRA_TOP_SW) THEN
!       Calculate the radiation absorbed in the extra layer
!       above the top of the rest of the model.
        DO L=1, NLIT
          TOP_ABSORPTION(LIST(L))=FLUX_NET(L, 0)                        &
     &      -FLUX_NET(L, N_LAYER-NLEVS)
        ENDDO
      ENDIF
!
!
! Extra direct SW flux reaching the surface due to the
! orography correction:

      IF (L_orog) THEN
         ALLOCATE(lg_f_orog(NLIT))

         lg_f_orog = FLUX_DIRECT(1:NLIT, N_LAYER) *                     &
     &           (lg_orog_corr - 1.0)/lg_orog_corr

         DO L=1, NLIT
           f_orog(col_list(L),row_list(L)) = lg_f_orog(L)
         ENDDO
      ENDIF
!
!
!     ASSIGNMENT OF DIAGNOSTICS:
!
!     Note: Purely diagnostic quantities allocated dynamically in
!     RAD_CTL2 are zeroed there and need to be filled only at lit
!     points.
!
!
!     OUTGOING SOLAR RADIATION AT TOA:
!
      IF (SW_diag%L_solar_out_toa) THEN
         DO L=1, NLIT
            SW_diag%solar_out_toa(col_list(l), row_list(l))             &
     &         =SOLAR_INCIDENT_NORM(L)/SEC_0(L)                         &
     &         -FLUX_NET(L, 0)
         ENDDO
      ENDIF
!
!
!     CLEAR-SKY OUTGOING SOLAR RADIATION AT TOA:
!
      IF (SW_diag%L_solar_out_clear) THEN
         DO L=1, NLIT
            SW_diag%solar_out_clear(col_list(l), row_list(l))           &
     &         =SOLAR_INCIDENT_NORM(L)/SEC_0(L)                         &
     &         -FLUX_NET_CLEAR(L, 0)
         ENDDO
      ENDIF
!
!
!     SURFACE FLUX BELOW 690NM.
!
      IF (L_FLUX_BELOW_690NM_SURF) THEN
        IF(L_CTILE)THEN
          DO L=1, NLIT
            IF (FLANDG(LIST(L)) <  1.0) THEN
              SW_diag%FlxSeaBelow690nmSurf(col_list(l), row_list(l))    &
     &          =FL_SEA_BELOW_690NM_SURF_G(L)                           &
     &         *(1.0E+00-ICE_FRACTION(LIST(L)))
            ENDIF
            FRACSOLID=FLANDG(LIST(L))                                   &
     &        +(1.-FLANDG(LIST(L)))*ICE_FRACTION(LIST(L))
            IF(FRACSOLID >  0.0)THEN
              SW_diag%FlxSolBelow690nmSurf(col_list(l), row_list(l))=   &
     &          (FLUX_BELOW_690NM_SURF_G(L)-                            &
     &          (1.-FLANDG(LIST(L)))*                                   &
     &          SW_diag%FlxSeaBelow690nmSurf(col_list(l), row_list(l))) &
     &          /FRACSOLID
            ENDIF
          ENDDO
        ELSE

          DO L=1, NLIT
            IF (LAND(LIST(L))) THEN
               FLUX_BELOW_690NM_SURF(LIST(L))                           &
     &            =FLUX_BELOW_690NM_SURF_G(L)
            ELSE
               FLUX_BELOW_690NM_SURF(LIST(L))                           &
     &            =FLUX_BELOW_690NM_SURF_G(L)                           &
     &            *(1.0E+00-ICE_FRACTION(LIST(L)))
            ENDIF
          ENDDO
        ENDIF

      ENDIF

!     Orography correction to direct SW flux:
!
      IF (SW_diag%L_orog_corr) THEN
        DO L=1, NLIT
          SW_diag%orog_corr(col_list(l),row_list(l))=lg_orog_corr(L)
        ENDDO
      ENDIF
!
!
!     COMPONENTS OF DOWNWARD FLUX AT THE SURFACE FOR MOSES II
!
      IF (L_MOSES_II .or. l_cable) THEN
        DO L=1, NLIT
           SURF_DOWN_SW(LIST(L),1) = SURF_VIS_DIR_G(L)
           SURF_DOWN_SW(LIST(L),2) = SURF_VIS_DIF_G(L)
           SURF_DOWN_SW(LIST(L),3) = SURF_NIR_DIR_G(L)
           SURF_DOWN_SW(LIST(L),4) = SURF_NIR_DIF_G(L)
        ENDDO
      ENDIF
!
!
!     DOWNWARD FLUX AT THE SURFACE:
!
      IF (SW_diag%L_surface_down_flux) THEN
         DO L=1, NLIT
            SW_diag%surface_down_flux(col_list(l), row_list(l))         &
     &        =SURFACE_DOWN_FLUX_G(L)
         ENDDO
      ENDIF
!
!
!     CLEAR-SKY DOWNWARD FLUX AT THE SURFACE:
!
      IF (SW_diag%L_surf_down_clr) THEN
         DO L=1, NLIT
            SW_diag%surf_down_clr(col_list(l), row_list(l))             &
     &        =SURF_DOWN_CLR_G(L)
         ENDDO
      ENDIF
!
!
!     CLEAR-SKY UPWARD FLUX AT THE SURFACE:
!
      IF (SW_diag%L_surf_up_clr) THEN
         DO L=1, NLIT
            SW_diag%surf_up_clr(col_list(l), row_list(l))               &
     &        =SURF_UP_CLR_G(L)
         ENDDO
      ENDIF
!
!
!     NET FLUX AT THE TROPOPAUSE:
!
      IF (SW_diag%L_net_flux_trop) THEN
         DO L=1, NLIT
            SW_diag%net_flux_trop(col_list(l), row_list(l))             &
     &         =FLUX_NET(L, N_LAYER+1-TRINDX(LIST(L)))
         ENDDO
      ENDIF
!
!
!     UPWARD FLUX AT THE TROPOPAUSE:
!
      IF (SW_diag%L_up_flux_trop) THEN
         DO L=1, NLIT
            SW_diag%up_flux_trop(col_list(l), row_list(l))              &
     &         =FLUX_UP(L, N_LAYER+1-TRINDX(LIST(L)))
         ENDDO
      ENDIF

!
!  DIRECT AND DIFFUSE DOWNWARD FLUX
!
      IF (SW_diag%L_flux_direct) THEN
         DO i=1, nlevs+1
            DO l=1, nlit
               SW_diag%flux_direct(col_list(l), row_list(l),i)          &
     &                 =flux_direct(l,nlevs+1-i)
            ENDDO
         ENDDO
      ENDIF
      
      IF (SW_diag%L_flux_diffuse) THEN
         DO i=1, nlevs+1
            DO l=1, nlit
               SW_diag%flux_diffuse(col_list(l), row_list(l),i)         &
     &                 =flux_diffuse(l,nlevs+1-i)
            ENDDO
         ENDDO
      ENDIF
      
!
!
!    CLOUD EXTINCTION DIAGNOSTICS
!
      IF (SW_diag%L_cloud_extinction) THEN
        DO I=1, NCLDS
         DO L=1, NLIT
            SW_diag%cloud_extinction                                    &
     &         (col_list(l), row_list(l), I)                            &
     &         =CLOUD_EXTINCTION_G(L, N_LAYER+1-I)
             SW_diag%cloud_weight_extinction                            &
     &         (col_list(l), row_list(l), I)                            &
     &         =CLOUD_WEIGHT_EXTINCTION_G(L, N_LAYER+1-I)
         ENDDO
        ENDDO
      ENDIF
!
      IF (SW_diag%L_ls_cloud_extinction) THEN
        DO I=1, NCLDS
         DO L=1, NLIT
            SW_diag%ls_cloud_extinction                                 &
     &         (col_list(l), row_list(l), I)                            &
     &         =LS_CLOUD_EXTINCTION_G(L, N_LAYER+1-I)
            SW_diag%ls_cloud_weight_extinction                          &
     &         (col_list(l), row_list(l), I)                            &
     &         =LS_CLOUD_WEIGHT_EXTINCTION_G(L, N_LAYER+1-I)
         ENDDO
        ENDDO
      ENDIF
!
      IF (SW_diag%L_cnv_cloud_extinction) THEN
        DO I=1, NCLDS
         DO L=1, NLIT
            SW_diag%cnv_cloud_extinction                                &
     &         (col_list(l), row_list(l), I)                            &
     &         =CNV_CLOUD_EXTINCTION_G(L, N_LAYER+1-I)
            SW_diag%cnv_cloud_weight_extinction                         &
     &         (col_list(l), row_list(l), I)                            &
     &         =CNV_CLOUD_WEIGHT_EXTINCTION_G(L, N_LAYER+1-I)
         ENDDO
        ENDDO
      ENDIF
!

!
!
!
!
!     FINAL PROCESSING OF OUTPUT FIELDS
!
      IF (L_orog) THEN
         SWOUT_TEMP=SWOUT(:,1)
         SWOUT(LIST(1:NLIT),1)=SWOUT(LIST(1:NLIT),1) - lg_f_orog

         DEALLOCATE(lg_f_orog)
         DEALLOCATE(lg_orog_corr)
      ENDIF


!     CONVERT THE FLUXES TO INCREMENTS.
      DO I=NLEVS, 1, -1
!
        IF (l_mixing_ratio) THEN

!       The layer_heat_capacity array has been calculated using the
!       specific heat of moist air and the true layer mass when
!       l_mixing_ratio is true.

          DO L=1, NLIT
            SWOUT(LIST(L), I+1)=(SWOUT(LIST(L), I+1)-SWOUT(LIST(L), I)) &
     &         *PTS/layer_heat_capacity(L, N_LAYER+1-I)
          ENDDO

          IF (SW_diag%L_clear_hr) THEN
            DO L=1, NLIT
               SW_diag%clear_hr(col_list(l), row_list(l), I)            &
     &            =(FLUX_NET_CLEAR(L, N_LAYER-I)                        &
     &            -FLUX_NET_CLEAR(L, N_LAYER+1-I))                      &
     &            /layer_heat_capacity(L, N_LAYER+1-I)
            ENDDO
          ENDIF
!
        ELSE

!       When l_mixing_ratio is false, layer_heat_capacity is based
!       on the specific heat of dry air and a mass derived from the 
!       hydrostatic approximation as used below. For bit-compatibility
!       reasons the original equations have been retained for this case.

          DO L=1, N_POINTS
            SWOUT(L, I+1)=(SWOUT(L, I+1)-SWOUT(L, I))                   &
     &         /((P_LAYER_BOUNDARIES(L,I-1) -                           &
     &         P_LAYER_BOUNDARIES(L,I))*CPBYG/PTS)
          ENDDO
!
          IF (SW_diag%L_clear_hr) THEN
            DO L=1, NLIT
               SW_diag%clear_hr(col_list(l), row_list(l), I)            &
     &            =(FLUX_NET_CLEAR(L, N_LAYER-I)                        &
     &            -FLUX_NET_CLEAR(L, N_LAYER+1-I))                      &
     &            /((P_LAYER_BOUNDARIES(LIST(L), I-1) -                 &
     &            P_LAYER_BOUNDARIES(LIST(L), I))*CPBYG)
            ENDDO
          ENDIF

        ENDIF
!
      ENDDO
!
      IF (L_orog) SWOUT(:,1)=SWOUT_TEMP
!
!
!     SEPARATE CONTRIBUTIONS OVER OPEN SEA.
!     SEA_FLUX_G IS NOT WEIGHTED BY THE FRACTION OF ICE.
! Weight open sea flux with open sea fraction over total sea
      IF(L_CTILE)THEN
        DO L=1, NLIT
          IF (FLANDG(LIST(L)) <  1.0) THEN
            SWSEA(LIST(L))=(1.0E+00-ICE_FRACTION(LIST(L)))              &
     &         *SEA_FLUX_G(L)
            SWOUT(LIST(L), 1)=SWOUT(LIST(L), 1)                         &
     &        -(1.0E+00-FLANDG(LIST(L)))*SWSEA(LIST(L))
           ELSE
            SWSEA(LIST(L))=0.0
           ENDIF
         ENDDO
       ENDIF
       IF(.NOT.L_CTILE)THEN
!DIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO L=1, NLIT
           IF (.NOT.LAND(LIST(L))) THEN
             SWSEA(LIST(L))=(1.0E+00-ICE_FRACTION(LIST(L)))             &
     &         *SEA_FLUX_G(L)
             SWOUT(LIST(L), 1)=SWOUT(LIST(L), 1)-SWSEA(LIST(L))
          ENDIF
        ENDDO
      ENDIF

!
!
!     TOTAL DOWNWARD FLUX OF PHOTOSYTHETICALLY ACTIVE RADIATION. ADD
!     THIS TO THE SWOUT ARRAY AS AN EXTRA 'LEVEL' TO ENABLE USE IN NON-
!     RADIATION TIMESTEPS.
!
! Set dirpar_out to zero if selected. Will only be defined if
! MOSESII is selecetd as well.
      IF (l_direct_par) dirpar_out = 0.0
      IF (L_FLUX_BELOW_690NM_SURF) THEN
        IF (L_MOSES_II .or. l_cable) THEN
          DO L=1, NLIT
            SWOUT(LIST(L),NLEVS+2) = SURF_VIS_DIR_G(L) +                &
     &                               SURF_VIS_DIF_G(L)
          ENDDO
          IF (l_direct_par) THEN
             DO L=1, NLIT
                dirpar_out(LIST(L)) = SURF_VIS_DIR_G(L)
             ENDDO
          ENDIF
        ELSE
          DO L=1, N_POINTS
             IF(LAND(L))THEN
               SWOUT(L, NLEVS+2)=FLUX_BELOW_690NM_SURF(L) /             &
     &           (1 - LAND_ALB(L))
              ELSE
               SWOUT(L, NLEVS+2)=FLUX_BELOW_690NM_SURF(L) /             &
     &           (1 - SICE_ALB(L))
              ENDIF

          ENDDO
        ENDIF
      ELSE
        DO L=1, N_POINTS
           SWOUT(L, NLEVS+2)=0.0
        ENDDO
      ENDIF
!
!
!     DIVIDE BY COSINE OF SOLAR ZENITH ANGLE TO PROVIDE VALUES FOR
!     UPPER ROUTINES. THIS APPLIES ONLY TO SWOUT. THE MACHINE TOLERANCE
!     IS ADDED TO MAINTAIN CONDITIONING.
      DO I=1, NLEVS+2
         DO L=1, N_POINTS
            SWOUT(L, I)=SWOUT(L, I)/(COSZIN(L)*LIT(L)+TINY(COSZIN))
         ENDDO
      ENDDO
      IF (l_direct_par) THEN
        DO l=1, n_points
          dirpar_out(l) = dirpar_out(l) /                               &
     &      (coszin(l)*lit(l) + tiny(coszin))
        END DO
      END IF
!
!
!     DIVIDE SURFACE DOWNWARD SW COMPONENTS BY COSINE OF SOLAR ZENITH
!     ANGLE FOR MOSES II.
      IF (L_MOSES_II .or. l_cable) THEN
        DO I=1, 4
          DO L=1, N_POINTS
            SURF_DOWN_SW(L, I) = SURF_DOWN_SW(L, I) /                   &
     &                          (COSZIN(L)*LIT(L)+TINY(COSZIN))
          ENDDO
        ENDDO
      ENDIF
!
!
      RETURN
      END SUBROUTINE R2_SWRAD
#endif
