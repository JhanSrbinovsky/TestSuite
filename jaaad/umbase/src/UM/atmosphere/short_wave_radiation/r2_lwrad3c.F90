#if defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Longwave Interface to the Edwards-Slingo Radiation Scheme.
!
! Purpose:
!   This routine prepares the call to the Edwards-Slingo radiation
!   scheme in the longwave.
!
! Method:
!   Principally, this routine transfers arrays into the correct formats.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_LWRAD3C(IERR                                        &
!                       Gaseous Mixing Ratios
     &   , H2O, CO2, O3                                                 &
     &   , CO2_DIM1, CO2_DIM2, CO2_3D, L_CO2_3D                         &
     &   , L_use_stochem_CH4, CH4_stochem                               &
! chemical greenhouse gas fields
     &   , ngrgas, grgas_field                                          &
     &   , N2O_MIX_RATIO, CH4_MIX_RATIO                                 &
     &   , CFC11_MIX_RATIO, CFC12_MIX_RATIO, CFC113_MIX_RATIO           &
     &   , HCFC22_MIX_RATIO, HFC125_MIX_RATIO, HFC134A_MIX_RATIO        &
!                       Thermodynamic Variables
     &   , TAC,  TSTAR,  TSTAR_SOLID, TSTAR_SEA, L_CTILE, PSTAR         &
     &   , P_LAYER_BOUNDARIES                                           &
     &   , P_LAYER_CENTRES                                              &
     &   , HEIGHT_THETA                                                 &
     &   , HEIGHT_RHO                                                   &
!                       Options for treating clouds
     &   , L_GLOBAL_CLOUD_TOP, GLOBAL_CLOUD_TOP                         &
     &   , L_INHOM_CLOUD, INHOM_CLOUD, DP_CORR_STRAT, DP_CORR_CONV      &
!                       Stratiform Cloud Fields
     &   , L_CLOUD_WATER_PARTITION, L_PC2                               &
     &   , LCA_AREA, LCA_BULK, LCCWC1, LCCWC2                           &
!                       Convective Cloud Fields
     &   , CCA, CCCWP, ccw, lcbase, CCB, CCT                            &
!                       Surface Fields
     &   , LAND, FLANDG, ICE_FRACTION                                   &
     &   , LYING_SNOW                                                   &
!                       Aerosol Fields
     &   , L_CLIMAT_AEROSOL, L_CLIM_AERO_HGT, L_HadGEM1_Clim_Aero       &
     &   , BL_DEPTH, N_LEVELS_BL                                        &
     &   , L_USE_CLEARRH, L_USE_DUST, DUST_DIM1, DUST_DIM2              &
     &   , DUST_1, DUST_2, DUST_3, DUST_4, DUST_5, DUST_6               &
     &   , L_USE_BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2, BIOGENIC       &
     &   , L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT                     &
     &   , SULP_DIM1,SULP_DIM2                                          &
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
     &   , LW_SPECTRUM                                                  &
!                       Algorithmic options
     &   , LW_CONTROL                                                   &
     &   , PTS, L_MOD_K_FLUX, L_WENYI, i_gather                         &
!                       Diagnostics
     &   , LW_diag, row_list, col_list                                  &
!                       Physical Dimensions
     &   , N_PROFILE, NLEVS, N_LAYER, NCLDS                             &
     &   , NWET, NOZONE, row_length, rows, NPD_FIELD                    &
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                           &
     &   , N_CCA_LEV                                                    &
!                       Output Fields
     &   , OLR, TOP_ABSORPTION, LWSEA, LWOUT                            &
           ! Variables needed to calculate layer masses
     &   , rho_r2, r_rho_levels, r_theta_levels                         &
     &   , q, qcl, qcf, qcf2, qrain, qgraup                             &
     &   , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
     &   )
!
!
!
!     Modules included.
!
      use DEC_SPEC
      use CONTROL_STRUC
      use lwrdiag_mod

      USE rad_switches_mod, ONLY:                                       &
          lrad_ccrad

      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED
#include "c_0_dg_c.h"
#include "c_g.h"
#include "c_r_cp.h"
!     INTERNAL DIMENSIONS OF THE CODE
#include "dimfix3a.h"
!     SPECTRAL REGIONS
#include "spcrg3a.h"
!     METHODS OF INTEGRATION
#include "angint3a.h"
!     METHODS OF SCATTERING
#include "sctmth3a.h"
!     NUMERICAL PRECISIONS
!     SOLVERS
#include "solver3a.h"
!     PHYSICAL CONSTANTS
#include "phycn03a.h"
!     UNIT NUMBERS FOR PRINTED OUTPUT
#include "stdio3a.h"
!     ERROR FLAGS
#include "error3a.h"
!    Components of clouds for two-stream radiation code.
#include "cldcmp3a.h"
!    Cloud types for two-stream radiation code.
#include "cldtyp3a.h"
!     Numbers for cloud overlap schemes
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
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
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
!
!     SPECTRAL DATA:
      TYPE (SPECTRUM) :: LW_SPECTRUM
!
!     Control data:
      TYPE (CONTROL_OPTION) ::  LW_CONTROL
!
!     GASEOUS MIXING RATIOS:
      REAL                                                              &
                !, INTENT(IN)
     &     H2O(NPD_FIELD, NWET)                                         &
!             MASS MIXING RATIO OF WATER
     &   , CO2                                                          &
!             MASS MIXING RATIO OF CO2
     &   , O3(NPD_FIELD, NOZONE)                                        &
!             MASS MIXING RATIOS OF OZONE
     &   , CH4_stochem(NPD_FIELD, NLEVS)                                &
!             Mass mixing ratio of CH4 from STOCHEM
     &   , N2O_MIX_RATIO                                                &
!             MASS MIXING RATIO OF NITROUS OXIDE
     &   , CH4_MIX_RATIO                                                &
!             MASS MIXING RATIO OF METHANE
     &   , CFC11_MIX_RATIO                                              &
!             MASS MIXING RATIO OF CFC11
     &   , CFC12_MIX_RATIO                                              &
!             MASS MIXING RATIO OF CFC12
     &   , CFC113_MIX_RATIO                                             &
!             MASS MIXING RATIO OF CFC113
     &   , HCFC22_MIX_RATIO                                             &
!             MASS MIXING RATIO OF HCFC22
     &   , HFC125_MIX_RATIO                                             &
!             MASS MIXING RATIO OF HFC125
     &   , HFC134A_MIX_RATIO
!             MASS MIXING RATIO OF HFC134A
      INTEGER, INTENT(IN) :: ngrgas
      REAL, INTENT(IN) :: grgas_field(npd_field, nlevs, ngrgas)
!
!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL                                                              &
                !, INTENT(IN)
     &     P_LAYER_BOUNDARIES(NPD_FIELD,0:NLEVS)                        &
!             PRESSURE AT BOUNDARIES OF LAYERS
     &   , P_LAYER_CENTRES(NPD_FIELD,0:NLEVS)                           &
!             PRESSURE AT CENTRES OF LAYERS
     &   , HEIGHT_THETA(NPD_FIELD, 0:NLEVS)                             &
     &   , HEIGHT_RHO(NPD_FIELD, NLEVS)                                 &
     &   , TAC(NPD_FIELD, NLEVS)
!             TEMPERATURES AT CENTRES OF LAYERS

      Integer, intent(in) :: lcbase(npd_field)
      Real,    intent(in) ::    ccw(npd_field,nwet)

      Logical, intent(in)::                                             &
     &     l_mcr_qcf2                                                   &
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
     &   , DP_CORR_CONV
!             Decorrelation pressure scale for convective cloud
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
!             LIQUID WATER CONTENTS (THESE ARE NOT USED DIRECTLY IN
!             THE RADIATION: THE TOTAL CONDENSED WATER CONTENT IS
!             REPARTITIONED USING FOCWWIL).
     &   , LCCWC2(NPD_FIELD, NCLDS+1/(NCLDS+1))                         &
!             ICE WATER CONTENTS (THESE ARE NOT USED DIRECTLY IN
!             THE RADIATION: THE TOTAL CONDENSED WATER CONTENT IS
!             REPARTITIONED USING FOCWWIL).
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
!             FRACTION OF GRID-BOX COVERED BY CONVECTIVE CLOUD

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
     &   , L_WENYI                                                      &
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
                              ! USE DIRECT RAD. EFFECT OF OCFF AER
     &   , L_USE_OCFF_INDIRECT                                          &
                                ! USE INDIRECT EFFECT OF OCFF AER
     &   , L_USE_SEASALT_INDIRECT                                       &
!             FLAG TO USE SEA-SALT FOR INDIRECT EFFECT
     &   , L_USE_SEASALT_DIRECT                                         &
!             FLAG TO USE SEA-SALT FOR DIRECT EFFECT
     &   , L_use_stochem_CH4
!             Flag to use stochem methane field
      INTEGER                                                           &
                !,INTENT (IN)
     &     SULP_DIM1,SULP_DIM2                                          &
!             DIMENSIONS FOR _SULPHATE ARRAYS, (P_FIELD,P_LEVELS or 1,1)
     &   , DUST_DIM1, DUST_DIM2                                         &
!        dimensions for mineral dust arrays (p_field,p_levels or 1,1)
     &   , BIOGENIC_DIM1, BIOGENIC_DIM2                                 &
!        dimensions for biogenic aerosol array (p_field,p_levels or 1,1)
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
!             MIXING RATIO OF BIOGENIC AEROSOL
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
!             Mass mixing ratio of fresh fossil-fuel organic carbon
     &   , AGED_OCFF(OCFF_DIM1, OCFF_DIM2)                              &
!             Mass mixing ratio of aged fossil-fuel organic carbon
     &   , CLOUD_OCFF(OCFF_DIM1, OCFF_DIM2)                             &
!             Mass mixing ratio of in-cloud fossil-fuel organic carbon
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
!     SURFACE FIELDS:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     LAND(NPD_FIELD)                                              &
!             LAND MASK (True if land fraction >0.5)
     &   , L_CTILE
!             COASTAL TILING SWITCH
      REAL                                                              &
                !, INTENT(IN)
     &     FLANDG(NPD_FIELD)
!            Land fraction in grid box
      REAL                                                              &
                !, INTENT(IN)
     &     PSTAR(NPD_FIELD)                                             &
!             SURFACE PRESSURES
     &   , TSTAR(NPD_FIELD)                                             &
!             SURFACE TEMPERATURES
     &   , TSTAR_SOLID(NPD_FIELD)                                       &
!             SOLID SURFACE TEMPERATURE
!             (I.E. AREAL MEAN OF LAND AND SEA-ICE)
     &   , TSTAR_SEA(NPD_FIELD)                                         &
!             OPEN SEA SURFACE TEMPERATURES
     &   , ICE_FRACTION(NPD_FIELD)                                      &
!             SEA ICE FRACTION OF SEA PORTION OF GRID BOX
     &   , LYING_SNOW(NPD_FIELD)                                        &
!             MASS LOADING OF LYING SNOW
     &   , BL_DEPTH(NPD_FIELD)
!             Depth of the boundary layer
!
      Real                                                              &
                !, Intent(IN)
     &     Ntot_land                                                    &
                               ! Number of droplets over land / m-3
     &   , Ntot_sea            ! Number of droplets over sea / m-3
!
!                       Level of tropopause
      INTEGER, INTENT(IN) ::                                            &
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

      Integer, Intent(IN) :: i_gather(NPD_FIELD)
!                              List of points where radiation is to be
!                              calculated
!
!
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

!     CALCULATED FLUXES:
      REAL                                                              &
                !, INTENT(OUT)
     &     OLR(NPD_FIELD)                                               &
!             NET OUTGOING RADIATION
     &   , TOP_ABSORPTION(NPD_FIELD)                                    &
!             Absorption in the extra radiative layer at the top
!             of the model
     &   , LWOUT(NPD_FIELD, NLEVS+1)                                    &
!             NET DOWNWARD FLUXES OR HEATING RATES
     &   , LWSEA(NPD_FIELD)
!             SEA-SURFACE COMPONENTS OF FLUX
!
!
!
!
!     Diagnostics:
      Type (StrLWDiag) :: LW_diag
!
      Integer, Intent(IN) :: row_list(npd_field)
!                              List of row indices of points
!                              to be treated
      Integer, Intent(IN) :: col_list(npd_field)
!                              List of column indices of points
!                              to be treated
!
!
!
!     LOCAL VARIABLES.
!
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
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
     &  , I_GAS_OVERLAP(LW_SPECTRUM%NPD_BAND)                           &
!             OVERLAPS IN EACH BAND
     &   , I_SCATTER_METHOD_BAND(LW_SPECTRUM%NPD_BAND)
!             Method of treating scattering in the calculation
!             of optical propeties in each band
!
!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL                                                              &
     &     D_MASS(NPD_PROFILE, NPD_LAYER)                               &
!             MASS THICKNESSES OF LAYERS
     &   , P(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             PRESSURE FIELD
     &   , T(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             TEMPERATURE FIELD
     &   , T_BDY(NPD_PROFILE, 0: NPD_LAYER)                             &
!             TEMPERATURE FIELD AT BOUNDARIES
     &   , GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                      &
     &                       ,LW_SPECTRUM%NPD_SPECIES)
!             MASS FRACTIONS OF GASES
!
      REAL :: layer_heat_capacity(npd_profile, npd_layer)
!             Specific heat capacity of layer * d_mass

!     Surface Fields:
      Logical :: Land_g(npd_profile)
!             Gathered land-surface mask
      INTEGER                                                           &
     &     I_SURFACE(NPD_PROFILE)
!             TYPE OF SURFACE AT THE FOOT OF EACH PROFILE
      Real :: ice_fraction_g(npd_profile)
!             Gathered ice fraction
      REAL                                                              &
     &     ALBEDO_FIELD_DIFF(NPD_PROFILE,LW_SPECTRUM%NPD_BAND)          &
!             DIFFUSE ALBEDOS
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE, LW_SPECTRUM%NPD_BAND)          &
!             DIRECT ALBEDOS
     &   , EMISSIVITY_FIELD(NPD_PROFILE, LW_SPECTRUM%NPD_BAND)          &
!             EMISSIVITIES
     &   , T_SOLID(NPD_PROFILE)                                         &
!             GATHERED TEMPERATURE OF SOLID SURFACE
     &   , T_SEA(NPD_PROFILE)                                           &
!             GATHERED OPEN SEA TEMPERATURE
     &   , FLANDG_G(NPD_PROFILE)                                        &
!             GATHERED LAND FRACTION
     &   , ALBEDO_SEA_DIFF(NPD_PROFILE, LW_SPECTRUM%NPD_BAND)           &
!             DIFFUSE ALBEDO OF OPEN SEA
     &   , ALBEDO_SEA_DIR(NPD_PROFILE, LW_SPECTRUM%NPD_BAND)            &
!             DIRECT ALBEDO OF OPEN SEA
     &   , T_SURFACE(NPD_PROFILE)
!             GATHERED TEMPERATURE OF SURFACE
!
!     CLOUDY PROPERTIES:
      INTEGER                                                           &
     &     N_CONDENSED                                                  &
!             NUMBER OF CONDENSED PHASES
     &   , TYPE_CONDENSED(NPD_CLOUD_COMPONENT)                          &
!             TYPES OF CONDENSED COMPONENTS
     &   , I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)                       &
!             PARAMETRIZATION SCHEMES FOR COMPONENTS
     &   , N_CLOUD_TOP_GLOBAL
!             INVERTED GLOBAL TOPMOST CLOUDY LAYER
      REAL                                                              &
     &     CONDENSED_PARAM_LIST(LW_SPECTRUM%NPD_CLOUD_PARAMETER         &
     &        , NPD_CLOUD_COMPONENT, LW_SPECTRUM%NPD_BAND)              &
!             PARAMETERS FOR CONDENSED PHASES
     &   , CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER                 &
     &        , NPD_CLOUD_COMPONENT)                                    &
!             CHARACTERISTIC DIMENSIONS OF CONDENSED SPECIES
     &   , CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                &
     &        , NPD_CLOUD_COMPONENT)                                    &
!             MASS FRACTIONS OF LIQUID WATER
     &   , W_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             CLOUD AMOUNTS
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             CLOUD AMOUNTS
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
     &        , LW_SPECTRUM%NPD_AEROSOL_SPECIES)
!             MIXING RATIOS OF AEROSOLS
!
!     COUPLING FIELDS:
      INTEGER                                                           &
     &     N_FRAC_SOL_POINT                                             &
!             NUMBER OF POINTS WITH FRACTIONAL ICE/LAND COVER
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)
!             INDICES OF POINTS WITH FRACTIONAL ICE/LAND COVER
!
!     FLUXES:
      REAL                                                              &
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DIRECT FLUX
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                 &
!             CLEAR-SKY DIRECT FLUX
     &   , FLUX_NET(NPD_PROFILE, 0: NPD_LAYER)                          &
!             DOWNWARD/NET FLUX
     &   , FLUX_NET_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                    &
!             CLEAR-SKY DOWNWARD/NET FLUX
     &   , FLUX_UP(NPD_PROFILE, 0: NPD_LAYER)                           &
!             UPWARD FLUX
     &   , FLUX_UP_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                     &
!             CLEAR-SKY UPWARD FLUX
     &   , SURFACE_DOWN_FLUX_G(NPD_PROFILE)                             &
!             DOWNWARD FLUX AT SURFACE
     &   , SURF_DOWN_CLR_G(NPD_PROFILE)                                 &
!             CLEAR-SKY DOWNWARD FLUX AT SURFACE
     &   , LWSEA_LOCAL(NPD_PROFILE)
!             Net downward flux into the sea at gathered points
!
!     LOCAL FIELDS FOR CLOUD DIAGNOSTICS:
      REAL                                                              &
     &     CLOUD_ABSORPTIVITY_G(NPD_PROFILE, NPD_LAYER)                 &
!             MEAN ABSORPTION COEFFICIENT IN CLOUDS WEIGHTED BY THE
!             CLOUD AMOUNT AND THE CLEAR-SKY FLUX
     &   , CLOUD_WEIGHT_ABSORPTIVITY_G(NPD_PROFILE, NPD_LAYER)          &
!             WEIGHTING FACTOR FOR ABSORPTION IN CLOUDS: THE PRODUCT
!             OF THE CLOUD AMOUNT AND THE CLEAR-SKY FLUX
     &   , LS_CLOUD_ABSORPTIVITY_G(NPD_PROFILE, NPD_LAYER)              &
!           MEAN ABSORPTION COEFFICIENT IN LAYER CLOUDS WEIGHTED BY THE
!           CLOUD AMOUNT AND THE CLEAR-SKY FLUX
     &   , LS_CLOUD_WEIGHT_ABSORPTIVITY_G(NPD_PROFILE, NPD_LAYER)       &
!           WEIGHTING FACTOR FOR ABSORPTION IN LAYER CLOUDS: THE PRODUCT
!           OF THE CLOUD AMOUNT AND THE CLEAR-SKY FLUX
     &   , CNV_CLOUD_ABSORPTIVITY_G(NPD_PROFILE, NPD_LAYER)             &
!           MEAN ABSORPTION COEFFICIENT IN CONV.CLOUDS WEIGHTED BY THE
!           CLOUD AMOUNT AND THE CLEAR-SKY FLUX
     &   , CNV_CLOUD_WEIGHT_ABSORPTIVITY_G(NPD_PROFILE, NPD_LAYER)
!           WEIGHTING FACTOR FOR ABSORPTION IN CONV.CLOUDS: THE PRODUCT
!           OF THE CLOUD AMOUNT AND THE CLEAR-SKY FLUX
!
!     Local Arrays to fill diagnostics:
      Real :: Total_cloud_cover_g(NPD_PROFILE)
!             Cloud fraction at gathered points
!
!     FIELDS REQUIRED FOR CALL TO RADIATION CODE BUT NOT USED
      INTEGER                                                           &
     &     N_ORDER_GAUSS                                                &
     &   , I_GAS
      REAL                                                              &
     &     SEC_0(NPD_PROFILE)                                           &
     &   , SOLAR_CONSTANT(NPD_PROFILE)
!
!
!     AUXILIARY VARIABLES:
      REAL                                                              &
     &     CPBYG                                                        &
!             SPECIFIC HEAT BY GRAVITY
     &   , WEIGHT_BAND(LW_SPECTRUM%NPD_BAND)                            &
!             WEIGHTING FACTORS FOR BANDS
     &   , NULLMMR
!             NULL MASS MIXING RATIO
      PARAMETER(CPBYG=CP/G)
      PARAMETER(NULLMMR=0.0E+00)
!
!     DUMMY FIELDS FOR RADIATION CODE
      LOGICAL                                                           &
     &     L_DUMMY
      REAL                                                              &
     &     DUMMY1D(1)                                                   &
     &    ,DUMMY2D(1,1)
!
!     LOCAL ARRAYS FOR AEROSOL OPTICAL DEPTH DIAGNOSTICS
      REAL                                                              &
     &     AOD_SULPHATE(NPD_PROFILE, LW_SPECTRUM%NPD_AOD_WAVEL),        &
     &     AOD_DUST(NPD_PROFILE, LW_SPECTRUM%NPD_AOD_WAVEL),            &
     &     AOD_SEASALT(NPD_PROFILE, LW_SPECTRUM%NPD_AOD_WAVEL),         &
     &     AOD_SOOT(NPD_PROFILE, LW_SPECTRUM%NPD_AOD_WAVEL),            &
     &     AOD_BIOMASS(NPD_PROFILE, LW_SPECTRUM%NPD_AOD_WAVEL),         &
     &     AOD_BIOGENIC(NPD_PROFILE, LW_SPECTRUM%NPD_AOD_WAVEL),        &
     &     AOD_OCFF(NPD_PROFILE, LW_SPECTRUM%NPD_AOD_WAVEL),            &
     &     AOD_DELTA(NPD_PROFILE, LW_SPECTRUM%NPD_AOD_WAVEL)
!
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     R2_SET_GAS_MIX_RATIO, R2_SET_THERMODYNAMIC                   &
     &   , R2_SET_AEROSOL_FIELD, R2_SET_CLOUD_FIELD                     &
     &   , R2_SET_CLOUD_PARAMETRIZATION                                 &
     &   , R2_SET_SURFACE_FIELD_LW                                      &
     &   , R2_COMPARE_PROC
!
!
!
!
!     INITIALIZE THE ERROR FLAG FOR THE RADIATION CODE.
      IERR=I_NORMAL
!     SET THE LOGICAL FLAG FOR DUMMY DIAGNOSTICS NOT AVAILABLE FROM
!     THE LOWER CODE IN THE LONG-WAVE TO .FALSE..
      L_DUMMY=.FALSE.

      IF (lrad_ccrad) THEN 
        CONDENSED_DIM_CHAR(:,:,:) = 0.0E+00
      END IF

!
!
!     COMPARE PROCESSES IN THE SPECTRAL FILE WITH THOSE ENABLED IN
!     THE CODE. ALSO WARN ABOUT COLLISION BETWEEN INTERACTIVE AEROSOLS
!     AND THE NWP AEROSOL CLIMATOLOGY.
! DEPENDS ON: r2_compare_proc
      CALL R2_COMPARE_PROC(IERR, LW_SPECTRUM%L_PRESENT                  &
     &   , LW_CONTROL%L_RAYLEIGH, LW_CONTROL%L_GAS                      &
     &   , LW_CONTROL%L_CONTINUUM, LW_CONTROL%L_DROP                    &
     &   , LW_CONTROL%L_AEROSOL, LW_CONTROL%L_AEROSOL_CCN               &
     &   , LW_CONTROL%L_ICE                                             &
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
     &   , LW_SPECTRUM%NPD_TYPE                                         &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!
!     SET THE MIXING RATIOS OF GASES.
! DEPENDS ON: r2_set_gas_mix_ratio
      CALL R2_SET_GAS_MIX_RATIO(IERR                                    &
     &   , N_PROFILE, NLEVS, N_LAYER, NWET, NOZONE                      &
     &   , I_GATHER                                                     &
     &   , LW_CONTROL%L_EXTRA_TOP                                       &
     &   , LW_SPECTRUM%N_ABSORB                                         &
     &   , LW_SPECTRUM%TYPE_ABSORB                                      &
     &   , LW_CONTROL%L_N2O                                             &
     &   , LW_CONTROL%L_CH4                                             &
     &   , LW_CONTROL%L_CFC11                                           &
     &   , LW_CONTROL%L_CFC12                                           &
     &   , .FALSE.                                                      &
     &   , LW_CONTROL%L_CFC113                                          &
     &   , LW_CONTROL%L_HCFC22                                          &
     &   , LW_CONTROL%L_HFC125                                          &
     &   , LW_CONTROL%L_HFC134A                                         &
     &   , H2O, CO2, O3, N2O_MIX_RATIO, CH4_MIX_RATIO                   &
     &   , CFC11_MIX_RATIO, CFC12_MIX_RATIO, NULLMMR                    &
     &   , CFC113_MIX_RATIO, HCFC22_MIX_RATIO, HFC125_MIX_RATIO         &
     &   , HFC134A_MIX_RATIO                                            &
     &   , GAS_MIX_RATIO                                                &
     &   , CO2_DIM1, CO2_DIM2, CO2_3D, L_CO2_3D                         &
     &   , CH4_stochem, L_use_stochem_CH4                               &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER                            &
     &   , LW_SPECTRUM%NPD_SPECIES                                      &
! chemical greenhouse gas fields
     &   , ngrgas, grgas_field                                          &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!
!     CALCULATE PRESSURES AND TEMPERATURES.
! DEPENDS ON: r2_set_thermodynamic
      CALL R2_SET_THERMODYNAMIC(N_PROFILE,NLEVS,N_LAYER,nwet,I_GATHER   &
     &   , LW_CONTROL%L_EXTRA_TOP,.TRUE.                                &
     &   , PSTAR, TSTAR, TSTAR_SOLID, TSTAR_SEA                         &
     &   , P_LAYER_BOUNDARIES                                           &
     &   , P_LAYER_CENTRES                                              &
     &   , HEIGHT_THETA                                                 &
     &   , HEIGHT_RHO                                                   &
     &   , TAC                                                          &
     &   , rho_r2, r_rho_levels, r_theta_levels                         &
     &   , q, qcl, qcf, qcf2, qrain, qgraup                             &
     &   , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
     &   , P(:,1:), T(:,1:), T_BDY, T_SURFACE, T_SOLID, T_SEA, D_MASS   &
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
     &      , N_PROFILE, NLEVS, N_LAYER                                 &
     &      , LW_SPECTRUM%N_AEROSOL                                     &
     &      , LW_SPECTRUM%TYPE_AEROSOL                                  &
     &      , I_GATHER                                                  &
     &      , LW_CONTROL%L_EXTRA_TOP                                    &
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
     &      , LAND, LYING_SNOW, PSTAR                                   &
     &      , P_LAYER_BOUNDARIES                                        &
     &      , TRINDX                                                    &
     &      , AEROSOL_MIX_RATIO                                         &
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER                         &
     &      , LW_SPECTRUM%NPD_AEROSOL_SPECIES                           &
     &      )
      ENDIF
!
!
!     ASSIGN THE PROPERTIES OF CLOUDS. A DUMMY ARRAY MUST BE PASSED
!     FOR THE MICROPHYSICAL DIAGNOSTICS SINCE THEY ARE NOT AVAILABLE
!     THROUGH STASH IN THE LONG-WAVE.
!
! DEPENDS ON: r2_set_cloud_parametrization
      CALL R2_SET_CLOUD_PARAMETRIZATION(IERR, LW_SPECTRUM%N_BAND        &
     &   , LW_CONTROL%I_ST_WATER, LW_CONTROL%I_CNV_WATER                &
     &   , LW_CONTROL%I_ST_ICE, LW_CONTROL%I_CNV_ICE                    &
     &   , LW_SPECTRUM%L_DROP_TYPE                                      &
     &   , LW_SPECTRUM%I_DROP_PARAMETRIZATION                           &
     &   , LW_SPECTRUM%DROP_PARAMETER_LIST                              &
     &   , LW_SPECTRUM%DROP_PARM_MIN_DIM, LW_SPECTRUM%DROP_PARM_MAX_DIM &
     &   , LW_SPECTRUM%L_ICE_TYPE                                       &
     &   , LW_SPECTRUM%I_ICE_PARAMETRIZATION                            &
     &   , LW_SPECTRUM%ICE_PARAMETER_LIST                               &
     &   , LW_SPECTRUM%ICE_PARM_MIN_DIM, LW_SPECTRUM%ICE_PARM_MAX_DIM   &
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST                      &
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM                         &
     &   , LW_SPECTRUM%NPD_BAND, LW_SPECTRUM%NPD_DROP_TYPE              &
     &   , LW_SPECTRUM%NPD_ICE_TYPE, LW_SPECTRUM%NPD_CLOUD_PARAMETER    &
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
!     The gathered land flag and land fraction
!     are required for cloud microphysics.
      Do l=1, n_profile
        land_g(l)=land(i_gather(l))
        FLANDG_G(l)=FLANDG(i_gather(l))
      Enddo
!
! DEPENDS ON: r2_set_cloud_field
      CALL R2_SET_CLOUD_FIELD(N_PROFILE, NLEVS, N_LAYER, NCLDS          &
     &   , I_GATHER                                                     &
     &   , P, T, D_MASS                                                 &
     &   , CCB, CCT, CCA, CCCWP, ccw, lcbase                            &
     &   , LCCWC1, LCCWC2, LCA_AREA, LCA_BULK                           &
     &   , L_PC2, LW_CONTROL%L_MICROPHYSICS, L_AEROSOL_CCN              &
     &   , SEA_SALT_FILM, SEA_SALT_JET                                  &
     &   , L_USE_SEASALT_INDIRECT, SALT_DIM_IND_A, SALT_DIM_IND_B       &
     &   , L_USE_BIOGENIC, BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2       &
     &   , SULP_DIM1, SULP_DIM2, ACCUM_SULPHATE, DISS_SULPHATE          &
     &   , AITKEN_SULPHATE, L_USE_BMASS_INDIRECT                        &
     &   , BMASS_DIM1, BMASS_DIM2, AGED_BMASS, CLOUD_BMASS              &
     &   , L_USE_OCFF_INDIRECT, OCFF_DIM1, OCFF_DIM2                    &
     &   , AGED_OCFF, CLOUD_OCFF                                        &
     &   , LYING_SNOW                                                   &
     &   , L_CLOUD_WATER_PARTITION, land_g, FLANDG_G                    &
     &    , LW_CONTROL%I_CLOUD_REPRESENTATION, I_CONDENSED_PARAM        &
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM                         &
     &   , N_CONDENSED, TYPE_CONDENSED                                  &
     &   , W_CLOUD, FRAC_CLOUD, LW_CONTROL%L_LOCAL_CNV_PARTITION        &
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                      &
!                       Microphysical Diagnostics are not available
!                       in this spectral region.
     &   , DUMMY2D, .FALSE., DUMMY2D, .FALSE.                           &
     &   , DUMMY2D, .FALSE., DUMMY2D, .FALSE.                           &
     &   , DUMMY2D, .FALSE.                                             &
     &   , DUMMY2D, .FALSE., DUMMY2D, .FALSE.                           &
     &   , DUMMY2D, .FALSE., DUMMY2D, .FALSE.                           &
     &   , DUMMY2D, .FALSE., DUMMY2D, .FALSE.                           &
     &   , col_list, row_list, row_length, rows                         &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER                            &
     &   , LW_SPECTRUM%NPD_AEROSOL_SPECIES                              &
     &   , N_CCA_LEV, Ntot_land, Ntot_sea                               &
     &   )
!
!
! DEPENDS ON: r2_set_surface_field_lw
      CALL R2_SET_SURFACE_FIELD_LW(                                     &
     &     N_PROFILE, i_gather                                          &
     &   , LW_SPECTRUM%N_BAND                                           &
     &   , I_SURFACE, LW_SPECTRUM%I_SPEC_SURFACE                        &
     &   , LW_SPECTRUM%L_SURFACE                                        &
     &   , EMISSIVITY_FIELD, ALBEDO_FIELD_DIR, ALBEDO_FIELD_DIFF        &
     &   , ALBEDO_SEA_DIR, ALBEDO_SEA_DIFF                              &
     &   , FLANDG                                                       &
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION             &
     &   , ice_fraction_g                                               &
     &   , npd_field, NPD_PROFILE                                       &
     &   , LW_SPECTRUM%NPD_BAND, LW_SPECTRUM%NPD_SURFACE                &
     &   )
!
!
!     CHECK THAT A VALID NUMBER HAS BEEN SUPPLIED FOR THE SOLVER.
      IF ( (LW_CONTROL%I_SOLVER /= IP_SOLVER_PENTADIAGONAL).AND.        &
     &     (LW_CONTROL%I_SOLVER /= IP_SOLVER_MIX_11).AND.               &
     &     (LW_CONTROL%I_SOLVER /= IP_SOLVER_MIX_APP_SCAT).AND.         &
     &     (LW_CONTROL%I_SOLVER /= IP_SOLVER_MIX_DIRECT).AND.           &
     &     (LW_CONTROL%I_SOLVER /= IP_SOLVER_MIX_DIRECT_HOGAN).AND.     &
     &     (LW_CONTROL%I_SOLVER /= IP_SOLVER_HOMOGEN_DIRECT).AND.       &
     &     (LW_CONTROL%I_SOLVER /= IP_SOLVER_TRIPLE).AND.               &
     &     (LW_CONTROL%I_SOLVER /= IP_SOLVER_TRIPLE_HOGAN).AND.         &
     &     (LW_CONTROL%I_SOLVER /= IP_SOLVER_TRIPLE_APP_SCAT)           &
     &   ) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: AN INVALID SOLVER HAS BEEN SELECTED '           &
     &      , 'IN THE LONGWAVE REGION.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     SET CLEAR-SKY CALCULATIONS.
      L_CLEAR=LW_diag%L_clear_olr.OR.                                   &
     &        LW_diag%L_surf_down_clr.OR.                               &
     &        LW_diag%L_clear_hr.OR.                                    &
     &       (LW_diag%L_cloud_absorptivity.AND.                         &
     &              LW_diag%L_cloud_weight_absorptivity)
!
      IF (L_CLEAR) THEN
!
!       SELECT A CLEAR-SKY SOLVER TO MATCH THE MAIN SOLVER.
        IF (LW_CONTROL%I_SOLVER == IP_SOLVER_PENTADIAGONAL) THEN
           I_SOLVER_CLEAR=IP_SOLVER_PENTADIAGONAL
        ELSE IF (LW_CONTROL%I_SOLVER == IP_SOLVER_MIX_11) THEN
           I_SOLVER_CLEAR=IP_SOLVER_PENTADIAGONAL
        ELSE IF (LW_CONTROL%I_SOLVER == IP_SOLVER_MIX_APP_SCAT) THEN
           I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
        ELSE IF (LW_CONTROL%I_SOLVER == IP_SOLVER_MIX_DIRECT) THEN
           I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
        ELSE IF (LW_CONTROL%I_SOLVER == IP_SOLVER_MIX_DIRECT_HOGAN) THEN
           I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
        ELSE IF (LW_CONTROL%I_SOLVER == IP_SOLVER_TRIPLE) THEN
           I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
        ELSE IF (LW_CONTROL%I_SOLVER == IP_SOLVER_TRIPLE_HOGAN) THEN
           I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
        ELSE IF (LW_CONTROL%I_SOLVER == IP_SOLVER_TRIPLE_APP_SCAT) THEN
           I_SOLVER_CLEAR=IP_SOLVER_HOMOGEN_DIRECT
        ENDIF
!
      ENDIF
!
!     CHECK THE CONSISTENCY OF CLOUD DIAGNOSTICS
!
      IF (LW_diag%L_cloud_absorptivity) THEN
         IF (.NOT.LW_diag%L_cloud_weight_absorptivity) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** Error: The cloud absorptivity'                      &
     &         , 'may be diagnosed only in conjunction'                 &
     &         , 'with the corresponding weights.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (LW_diag%L_ls_cloud_absorptivity) THEN
         IF (.NOT.LW_diag%L_ls_cloud_weight_absorptivity) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** Error: The layer cloud absorptivity'                &
     &         , 'may be diagnosed only in conjunction'                 &
     &         , 'with the corresponding weights.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (LW_diag%L_cnv_cloud_absorptivity) THEN
         IF (.NOT.LW_diag%L_cnv_cloud_weight_absorptivity) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** Error: The conv. cloud absorptivity'                &
     &         , 'may be diagnosed only in conjunction'                 &
     &         , 'with the corresponding weights.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
!
!     SET PROPERTIES OF INDIVIDUAL BANDS.
      DO I=1, LW_SPECTRUM%N_BAND
         WEIGHT_BAND(I)=1.0E+00
         I_GAS_OVERLAP(I)=LW_CONTROL%I_GAS_OVERLAP
      ENDDO
!
!     Apply a consistent treatment of scattering across all bands
!     (This could be generalized if required.)
      DO I=1, LW_SPECTRUM%N_BAND
         I_SCATTER_METHOD_BAND(I)=LW_CONTROL%I_SCATTER_METHOD
      ENDDO
!
!     INVERT THE TOPMOST CLOUDY LAYER IF USING A GLOBAL VALUE.
      IF (L_GLOBAL_CLOUD_TOP) THEN
         N_CLOUD_TOP_GLOBAL=N_LAYER+1-GLOBAL_CLOUD_TOP
      ENDIF
!
!
!
! DEPENDS ON: flux_calc
      CALL FLUX_CALC(IERR                                               &
!                       Logical Flags for Processes
     &   , L_RAYLEIGH, L_AEROSOL, L_GAS, L_CONTINUUM                    &
     &   , LW_CONTROL%L_CLOUD, L_DROP, L_ICE, L_PC2                     &
!                       Angular Integration
     &   , LW_CONTROL%I_ANGULAR_INTEGRATION, LW_CONTROL%I_2STREAM       &
     &   , LW_CONTROL%L_2_STREAM_CORRECT                                &
     &   , LW_CONTROL%L_RESCALE, N_ORDER_GAUSS                          &
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_BAND                                        &
!                       Options for treating clouds
     &   , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL                       &
     &   , L_INHOM_CLOUD, INHOM_CLOUD                                   &
!                       Options for Solver
     &   , LW_CONTROL%I_SOLVER                                          &
!                       General Spectral Properties
     &   , LW_SPECTRUM%N_BAND, 1, LW_SPECTRUM%N_BAND, WEIGHT_BAND       &
!                       General Atmospheric Properties
     &   , N_PROFILE, N_LAYER                                           &
     &   , LW_CONTROL%L_LAYER, LW_CONTROL%L_CLOUD_LAYER                 &
     &   , P, T, T_SURFACE, T_SOLID, T_SEA, T_BDY, D_MASS               &
!                       Spectral Region
     &   , LW_CONTROL%ISOLIR                                            &
!                       Solar Fields
     &   , SEC_0, SOLAR_CONSTANT, LW_SPECTRUM%SOLAR_FLUX_BAND           &
     &   , LW_SPECTRUM%RAYLEIGH_COEFFICIENT                             &
!                       Infra-red Fields
     &   , LW_SPECTRUM%N_DEG_FIT                                        &
     &   , LW_SPECTRUM%THERMAL_COEFFICIENT                              &
     &   , LW_SPECTRUM%T_REF_PLANCK, LW_CONTROL%L_IR_SOURCE_QUAD        &
!                       Gaseous Absorption
     &   , LW_SPECTRUM%N_ABSORB, I_GAS_OVERLAP, I_GAS, GAS_MIX_RATIO    &
     &   , LW_SPECTRUM%N_BAND_ABSORB, LW_SPECTRUM%INDEX_ABSORB          &
     &   , LW_SPECTRUM%I_BAND_ESFT                                      &
     &   , LW_SPECTRUM%W_ESFT, LW_SPECTRUM%K_ESFT                       &
     &   , LW_SPECTRUM%I_SCALE_ESFT, LW_SPECTRUM%I_SCALE_FNC            &
     &   , L_WENYI, LW_SPECTRUM%SCALE_VECTOR                            &
     &   , LW_SPECTRUM%P_REFERENCE, LW_SPECTRUM%T_REFERENCE             &
     &   , L_MOD_K_FLUX                                                 &
!                       Doppler Broadening
     &   , LW_SPECTRUM%L_DOPPLER_PRESENT                                &
     &   , LW_SPECTRUM%DOPPLER_CORRECTION                               &
!                       Surface Fields
     &   , LW_SPECTRUM%L_SURFACE, I_SURFACE                             &
     &   , LW_SPECTRUM%I_SPEC_SURFACE                                   &
     &   , LW_SPECTRUM%SURFACE_ALBEDO                                   &
     &   , ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR                          &
     &   , LW_SPECTRUM%N_DIR_ALBEDO_FIT                                 &
     &   , LW_SPECTRUM%DIRECT_ALBEDO_PARM                               &
     &   , LW_SPECTRUM%EMISSIVITY_GROUND                                &
     &   , EMISSIVITY_FIELD                                             &
!                       Continuum Absorption
     &   , LW_SPECTRUM%N_BAND_CONTINUUM                                 &
     &   , LW_SPECTRUM%INDEX_CONTINUUM, LW_SPECTRUM%INDEX_WATER         &
     &   , LW_SPECTRUM%K_CONTINUUM                                      &
     &   , LW_SPECTRUM%I_SCALE_FNC_CONT                                 &
     &   , LW_SPECTRUM%SCALE_CONTINUUM                                  &
     &   , LW_SPECTRUM%P_REF_CONTINUUM                                  &
     &   , LW_SPECTRUM%T_REF_CONTINUUM                                  &
!                       Properties of Aerosols
     &   , LW_SPECTRUM%N_AEROSOL                                        &
     &   , AEROSOL_MIX_RATIO                                            &
     &   , LW_SPECTRUM%AEROSOL_ABSORPTION                               &
     &   , LW_SPECTRUM%AEROSOL_SCATTERING                               &
     &   , LW_SPECTRUM%AEROSOL_PHASE_FNC(:,1,:,:)                       &
     &   , LW_SPECTRUM%I_AEROSOL_PARAMETRIZATION                        &
     &   , LW_SPECTRUM%NHUMIDITY, LW_SPECTRUM%HUMIDITIES                &
     &   , LW_SPECTRUM%TYPE_AEROSOL, L_USE_CLEARRH                      &
!                       Optical depth of aerosols
     &   , LW_SPECTRUM%N_AOD_WAVEL                                      &
     &   , LW_diag%L_aod_sulphate, AOD_SULPHATE                         &
     &   , LW_diag%L_aod_dust,     AOD_DUST                             &
     &   , LW_diag%L_aod_seasalt,  AOD_SEASALT                          &
     &   , LW_diag%L_aod_soot,     AOD_SOOT                             &
     &   , LW_diag%L_aod_biomass,  AOD_BIOMASS                          &
     &   , LW_diag%L_aod_biogenic, AOD_BIOGENIC                         &
     &   , LW_diag%L_aod_ocff,     AOD_OCFF                             &
     &   , LW_diag%L_aod_delta,    AOD_DELTA                            &
     &   , LW_SPECTRUM%AOD_ABSORPTION, LW_SPECTRUM%AOD_SCATTERING       &
     &   , LW_SPECTRUM%I_AOD_TYPE                                       &
!                       Properties of Clouds
     &   , N_CONDENSED, TYPE_CONDENSED                                  &
     &   , LW_CONTROL%I_CLOUD, LW_CONTROL%I_CLOUD_REPRESENTATION        &
     &   , W_CLOUD, FRAC_CLOUD, TOT_CLOUD_COVER                         &
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                      &
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST                      &
     &   , DP_CORR_STRAT, DP_CORR_CONV                                  &
!                       Fluxes Calculated
     &   , FLUX_DIRECT, DUMMY2D, FLUX_NET, FLUX_UP                      &
     &   , DUMMY2D, DUMMY2D, DUMMY2D                                    &
     &   , L_DUMMY, L_DUMMY, L_DUMMY, L_DUMMY                           &
!                       Options for Clear-sky Fluxes
     &   , L_CLEAR, I_SOLVER_CLEAR                                      &
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_NET_CLEAR, FLUX_UP_CLEAR             &
!                       Arrays specific to the UM
!                       Arrays for Coupling
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION_G           &
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR, FLANDG_G, LWSEA_LOCAL       &
!                       Arrays for diagnostics specific to the UM
     &   , L_DUMMY, DUMMY1D, DUMMY1D, DUMMY1D, DUMMY1D                  &
     &   , L_DUMMY, L_CTILE, DUMMY1D, DUMMY1D, DUMMY1D, DUMMY1D         &
     &   , LW_diag%L_surface_down_flux, surface_down_flux_g             &
     &   , LW_diag%L_surf_down_clr, surf_down_clr_g                     &
     &   , L_DUMMY, DUMMY1D                                             &
     &   , L_DUMMY, DUMMY2D, DUMMY2D                                    &
     &   , L_DUMMY, DUMMY2D, DUMMY2D                                    &
     &   , L_DUMMY, DUMMY2D, DUMMY2D                                    &
     &   , LW_diag%L_cloud_absorptivity, CLOUD_absorptivity_G           &
     &   , CLOUD_WEIGHT_absorptivity_G                                  &
     &   , LW_diag%L_ls_cloud_absorptivity, LS_CLOUD_absorptivity_G     &
     &   , LS_CLOUD_WEIGHT_absorptivity_G                               &
     &   , LW_diag%L_cnv_cloud_absorptivity, CNV_CLOUD_absorptivity_G   &
     &   , CNV_CLOUD_WEIGHT_absorptivity_G                              &
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                           &
     &   , LW_SPECTRUM%NPD_BAND, LW_SPECTRUM%NPD_SPECIES                &
     &   , LW_SPECTRUM%NPD_ESFT_TERM, LW_SPECTRUM%NPD_SCALE_FNC         &
     &   , LW_SPECTRUM%NPD_SCALE_VARIABLE, LW_SPECTRUM%NPD_CONTINUUM    &
     &   , LW_SPECTRUM%NPD_AEROSOL_SPECIES, LW_SPECTRUM%NPD_HUMIDITIES  &
     &   , LW_SPECTRUM%NPD_CLOUD_PARAMETER                              &
     &   , LW_SPECTRUM%NPD_THERMAL_COEFF                                &
     &   , LW_SPECTRUM%NPD_SURFACE, LW_SPECTRUM%NPD_ALBEDO_PARM         &
     &   , LW_SPECTRUM%NPD_AOD_WAVEL                                    &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!
!
!     ASSIGNMENT OF DIAGNOSTICS:
!
!
!     OLR:
!
      DO L=1, N_PROFILE
         OLR(i_gather(l))=-FLUX_NET(L, 0)
      ENDDO
      IF (LW_diag%L_clear_olr) THEN
        Do l=1, n_profile
          LW_diag%clear_olr(col_list(l), row_list(l))                   &
     &      =-FLUX_NET_CLEAR(l, 0)
        Enddo
      ENDIF
!
!
!     TOTAL CLOUD COVER:
!
      IF (LW_diag%L_total_cloud_cover) THEN
        IF ( (LW_CONTROL%I_CLOUD == IP_CLOUD_PART_CORR).OR.             &
     &       (LW_CONTROL%I_CLOUD == IP_CLOUD_PART_CORR_CNV) ) THEN
          Do l=1, n_profile
            LW_diag%total_cloud_cover(col_list(l), row_list(l))         &
     &        =TOT_CLOUD_COVER(l)
          Enddo
        ELSE
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_CALC_TOTAL_CLOUD_COVER(N_PROFILE, N_LAYER, NCLDS      &
     &        , LW_CONTROL%I_CLOUD, W_CLOUD, total_cloud_cover_g        &
     &        , NPD_PROFILE, NPD_LAYER                                  &
     &        )
          Do l=1, n_profile
            LW_diag%total_cloud_cover(col_list(l), row_list(l))         &
     &        =total_cloud_cover_g(l)
          Enddo
        ENDIF
      ENDIF
!
!
!     NET FLUX AT THE TROPOPAUSE:
!
      IF (LW_diag%L_net_flux_trop) THEN
        Do l=1, n_profile
          LW_diag%net_flux_trop(col_list(l), row_list(l))               &
     &      =FLUX_NET(L, N_LAYER+1-TRINDX(i_gather(l)))
        Enddo
      ENDIF
!
!
!     DOWNWARD FLUX AT THE TROPOPAUSE:
!
      IF (LW_diag%L_down_flux_trop) THEN
        Do l=1, n_profile
          LW_diag%down_flux_trop(col_list(l), row_list(l))              &
     &      =FLUX_NET(L, N_LAYER+1-TRINDX(i_gather(l)))                 &
     &      +FLUX_UP(L, N_LAYER+1-TRINDX(i_gather(l)))
        Enddo
      ENDIF
!
!
!     Surface Fluxes:
      IF (LW_diag%L_surface_down_flux) THEN
        Do l=1, n_profile
          LW_diag%surface_down_flux(col_list(l), row_list(l))           &
     &      =surface_down_flux_g(l)
        Enddo
      Endif
      IF (LW_diag%L_surf_down_clr) THEN
        Do l=1, n_profile
          LW_diag%surf_down_clr(col_list(l), row_list(l))               &
     &      =surf_down_clr_g(l)
        Enddo
      Endif
!
!
!    CLOUD ABSORPTIVITY DIAGNOSTICS
!
      IF (LW_diag%L_cloud_absorptivity) THEN
        DO I=1, NCLDS
         DO L=1, n_profile
            LW_diag%cloud_absorptivity(col_list(l), row_list(l), I)     &
     &         =CLOUD_absorptivity_G(L, N_LAYER+1-I)
             LW_diag%cloud_weight_absorptivity(col_list(l),             &
     &         row_list(l), I)                                          &
     &         =CLOUD_WEIGHT_absorptivity_G(L, N_LAYER+1-I)
         ENDDO
        ENDDO
      ENDIF
!
      IF (LW_diag%L_ls_cloud_absorptivity) THEN
        DO I=1, NCLDS
         DO L=1, n_profile
            LW_diag%ls_cloud_absorptivity(col_list(l), row_list(l), I)  &
     &         =LS_CLOUD_absorptivity_G(L, N_LAYER+1-I)
            LW_diag%ls_cloud_weight_absorptivity(col_list(l),           &
     &         row_list(l), I)                                          &
     &         =LS_CLOUD_WEIGHT_absorptivity_G(L, N_LAYER+1-I)
         ENDDO
        ENDDO
      ENDIF
!
      IF (LW_diag%L_cnv_cloud_absorptivity) THEN
        DO I=1, NCLDS
         DO L=1, n_profile
            LW_diag%cnv_cloud_absorptivity(col_list(l), row_list(l), I) &
     &         =CNV_CLOUD_absorptivity_G(L, N_LAYER+1-I)
            LW_diag%cnv_cloud_weight_absorptivity(col_list(l),          &
     &         row_list(l), I)                                          &
     &         =CNV_CLOUD_WEIGHT_absorptivity_G(L, N_LAYER+1-I)
         ENDDO
        ENDDO
      ENDIF
!
!    TOTAL CLOUD FRACTION ON MODEL LEVELS
!
      IF (LW_diag%L_total_cloud_on_levels) THEN
        DO I=1, NCLDS
         DO L=1, n_profile
            LW_diag%total_cloud_on_levels(col_list(l), row_list(l), I)  &
     &         =W_CLOUD(L, N_LAYER+1-I)
         ENDDO
        ENDDO
      ENDIF
!
! ######################################################
! CLOUD WATER MIXING RATIOS
! ######################################################
!
!     LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (LIQUID)
!     ==============================================================
      If (LW_diag%L_ls_qcl_rad) Then
        Do I=1, NLEVS
          Do L=1, n_profile
            LW_diag%ls_qcl_rad(col_list(L),row_list(L),I)               &
     &         = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_ST_WATER)   &
     &         * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SW)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (ICE)
!     ==============================================================
      If (LW_diag%L_ls_qcf_rad) Then
        Do I=1, NLEVS
          Do L=1, n_profile
            LW_diag%ls_qcf_rad(col_list(L),row_list(L),I)               &
     &         = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_ST_ICE)     &
     &         * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SI)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (LIQUID)
!     ==============================================================
      If (LW_diag%L_cc_qcl_rad) Then
        Do I=1, NLEVS
          Do L=1, n_profile
            LW_diag%cc_qcl_rad(col_list(L),row_list(L),I)               &
     &         = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_CNV_WATER)  &
     &         * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CW)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (ICE)
!     ==============================================================

      If (LW_diag%L_cc_qcf_rad) Then
        Do I=1, NLEVS
          Do L=1, n_profile
            LW_diag%cc_qcf_rad(col_list(L),row_list(L),I)               &
     &         = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_CNV_ICE)    &
     &         * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CI)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
! ######################################################
! CLOUD FRACTIONS
! ######################################################
!
!     LARGE-SCALE cloud GRIDBOX FRACTION seen by radiation. (LIQUID)
!     ==============================================================
      If (LW_diag%L_ls_cl_rad) Then
        Do I=1, NLEVS
          Do L=1, n_profile
            LW_diag%ls_cl_rad(col_list(L),row_list(L),I)                &
     &         = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SW)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     LARGE-SCALE cloud GRIDBOX fraction seen by radiation. (ICE)
!     ==============================================================
      If (LW_diag%L_ls_cf_rad) Then
        Do I=1, NLEVS
          Do L=1, n_profile
            LW_diag%ls_cf_rad(col_list(L),row_list(L),I)                &
     &         = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SI)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     CONVECTIVE cloud GRIDBOX fraction seen by radiation. (LIQUID)
!     ==============================================================
      If (LW_diag%L_cc_cl_rad) Then
        Do I=1, NLEVS
          Do L=1, n_profile
            LW_diag%cc_cl_rad(col_list(L),row_list(L),I)                &
     &         = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CW)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     CONVECTIVE cloud GRIDBOX FRACTION seen by radiation. (ICE)
!     ==============================================================
      If (LW_diag%L_cc_cf_rad) Then
        Do I=1, NLEVS
          Do L=1, n_profile
            LW_diag%cc_cf_rad(col_list(L),row_list(L), I)               &
     &         = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CI)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif      
!
!    AEROSOL OPTICAL DEPTH
!
      IF (LW_diag%L_aod_sulphate) THEN
        DO I=1, LW_SPECTRUM%N_AOD_WAVEL
          DO L=1, n_profile
            LW_diag%aod_sulphate(col_list(L), row_list(L), I)           &
     &       = AOD_SULPHATE(L, I)
          ENDDO
        ENDDO
      ENDIF
      IF (LW_diag%L_aod_dust) THEN
        DO I=1, LW_SPECTRUM%N_AOD_WAVEL
          DO L=1, n_profile
            LW_diag%aod_dust(col_list(L), row_list(L), I)               &
     &       = AOD_DUST(L, I)
          ENDDO
        ENDDO
      ENDIF
      IF (LW_diag%L_aod_seasalt) THEN
        DO I=1, LW_SPECTRUM%N_AOD_WAVEL
          DO L=1, n_profile
            LW_diag%aod_seasalt(col_list(L), row_list(L), I)            &
     &       = AOD_SEASALT(L, I)
          ENDDO
        ENDDO
      ENDIF
      IF (LW_diag%L_aod_soot) THEN
        DO I=1, LW_SPECTRUM%N_AOD_WAVEL
          DO L=1, n_profile
            LW_diag%aod_soot(col_list(L), row_list(L), I)               &
     &       = AOD_SOOT(L, I)
          ENDDO
        ENDDO
      ENDIF
      IF (LW_diag%L_aod_biomass) THEN
        DO I=1, LW_SPECTRUM%N_AOD_WAVEL
          DO L=1, n_profile
            LW_diag%aod_biomass(col_list(L), row_list(L), I)            &
     &       = AOD_BIOMASS(L, I)
          ENDDO
        ENDDO
      ENDIF
      IF (LW_diag%L_aod_biogenic) THEN
        DO I=1, LW_SPECTRUM%N_AOD_WAVEL
          DO L=1, n_profile
            LW_diag%aod_biogenic(col_list(L), row_list(L), I)           &
     &       = AOD_BIOGENIC(L, I)
          ENDDO
        ENDDO
      ENDIF
      IF (LW_diag%L_aod_ocff) THEN
        DO I=1, LW_SPECTRUM%N_AOD_WAVEL
          DO L=1, n_profile
            LW_diag%aod_ocff(col_list(L), row_list(L), I)               &
     &       = AOD_OCFF(L, I)
          ENDDO
        ENDDO
      ENDIF
      IF (LW_diag%L_aod_delta) THEN
        DO I=1, LW_SPECTRUM%N_AOD_WAVEL
          DO L=1, n_profile
            LW_diag%aod_delta(col_list(L), row_list(L), I)              &
     &       = AOD_DELTA(L, I)
          ENDDO
        ENDDO
      ENDIF
!
!
!
!     OUTPUT ARRAYS:
!
!     CONVERT THE FLUXES TO INCREMENTS IN THE HEATING RATE EXCEPT AT
!     THE SURFACE: THERE, THE NET DOWNWARD FLUX IS ASSIGNED TO LWOUT.
      DO K=NLEVS, 1, -1

        IF (l_mixing_ratio) THEN
!       The layer_heat_capacity array has been calculated using the
!       specific heat of moist air and the true layer mass when
!       l_mixing_ratio is true.

          DO L=1, N_PROFILE
            LWOUT(i_gather(l), K+1)=(FLUX_NET(L, N_LAYER-K)             &
     &        -FLUX_NET(L, N_LAYER+1-K))                                &
     &        *PTS/layer_heat_capacity(L, N_LAYER+1-K)
          ENDDO
          IF (LW_diag%L_clear_hr) THEN
!           THE FACTOR OF PTS IS INCLUDED HERE TO YIELD A RATE FROM AN
!           INCREMENT.
            Do l=1, n_profile
              LW_diag%clear_hr(col_list(l), row_list(l), K)             &
     &          =(FLUX_NET_CLEAR(L, N_LAYER-K)                          &
     &          -FLUX_NET_CLEAR(L, N_LAYER+1-K))                        &
     &          /layer_heat_capacity(L, N_LAYER+1-K)
            Enddo
          ENDIF

        ELSE
!       When l_mixing_ratio is false, layer_heat_capacity is based
!       on the specific heat of dry air and a mass derived from the 
!       hydrostatic approximation as used below. For bit-compatibility
!       reasons the original equations have been retained for this case.

          DO L=1, N_PROFILE
            LWOUT(i_gather(l), K+1)=(FLUX_NET(L, N_LAYER-K)             &
     &        -FLUX_NET(L, N_LAYER+1-K))/                               &
     &        ((P_LAYER_BOUNDARIES(i_gather(l), K-1)                    &
     &        - P_LAYER_BOUNDARIES(i_gather(l), K))*CPBYG/PTS)
          ENDDO
          IF (LW_diag%L_clear_hr) THEN
!           THE FACTOR OF PTS IS INCLUDED HERE TO YIELD A RATE FROM AN
!           INCREMENT.
            Do l=1, n_profile
              LW_diag%clear_hr(col_list(l), row_list(l), K)             &
     &          =(FLUX_NET_CLEAR(L, N_LAYER-K)                          &
     &          -FLUX_NET_CLEAR(L, N_LAYER+1-K))/                       &
     &          ((P_LAYER_BOUNDARIES(i_gather(l), K-1) -                &
     &          P_LAYER_BOUNDARIES(i_gather(l), K))*CPBYG)
            Enddo
          ENDIF
        ENDIF
      ENDDO
!
      IF (LW_CONTROL%L_EXTRA_TOP) THEN
!       Calculate the radiation absorbed in the extra layer
!       above the top of the rest of the model.
        DO L=1, N_PROFILE
          TOP_ABSORPTION(i_gather(l))=FLUX_NET(L, 0)                    &
     &      -FLUX_NET(L, N_LAYER-NLEVS)
        ENDDO
      ENDIF

      DO L=1, N_PROFILE
         LWOUT(i_gather(l), 1)=FLUX_NET(L, N_LAYER)
      ENDDO
!
!     SEPARATE THE CONTRIBUTIONS OVER OPEN SEA AND SEA-ICE.
!     LWSEA MUST BE WEIGHTED WITH THE FRACTION OF OPEN SEA.
      DO L=1, N_PROFILE
         LWSEA(I_GATHER(L))=LWSEA_LOCAL(L)
         IF (FLANDG(I_GATHER(L)) == 1.0.OR.                             &
     &                         ICE_FRACTION(I_GATHER(L)) == 1.0) THEN
            LWSEA(I_GATHER(L))=0.0
         ELSE IF (FLANDG(I_GATHER(L)) <  EPSILON(LWOUT).AND.            &
     &      ICE_FRACTION(I_GATHER(L)) <  EPSILON(LWOUT)) THEN
            LWSEA(I_GATHER(L))=LWOUT(I_GATHER(L), 1)
            LWOUT(I_GATHER(L), 1)=0.0
         ELSE
!           LWSEA MUST BE SCALED BY THE FRACTION OF OPEN SEA TO
!           TOTAL SEA FOR CONSISTENCY WITH UPPER LEVELS IN THE MODEL.
            LWSEA(i_gather(l))=(1.0E+00-ICE_FRACTION(i_gather(l)))      &
     &        *LWSEA_LOCAL(l)
            LWOUT(I_GATHER(L), 1)=LWOUT(I_GATHER(L), 1)-                &
     &             (1.0E+00-FLANDG(I_GATHER(L)))*LWSEA(I_GATHER(L))
         ENDIF
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE R2_LWRAD3C
#endif
