
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Perform a 1-timestep integration of the Atmosphere Model
!
! Subroutine Interface:
      SUBROUTINE Atm_Step(                                              &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
     &  A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
     &  A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
     &  A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGLNDM Constants for physics routines
     &  land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Start arg_atm_fields.h

! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!
! Current Code Owner: A. Treshansky
!


      ! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP, &
      ! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV, &
      ! 1.2: Data variables stored in secondary space.
       P, & 
      ! Pressure on theta levels
       P_THETA_LEVELS, &
      ! Exner pressure on theta levels
       EXNER_THETA_LEVELS, &
      ! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN, &
      ! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF, &
      ! 1.5: Radiation Increments
       SW_INCS, LW_INCS, &
! PAR radiation increment
       DIRPAR, &
      ! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP, &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3, &
!  tropopause-based ozone
       TPPSOZONE, &
      ! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK, &       
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6, &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3, &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD, &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, &
       CO2, CH4_STOCH, O3_STOCH, &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5, &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10, &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15, &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20, &
! CABLE
!       TSOIL_TILE, SMCL_TILE, STHU_TILE, STHF_TILE, SNOW_DEPTH3L,       &
       TSOIL_TILE, SMCL_TILE,  STHF_TILE, SNOW_DEPTH3L,                 &
       SNOW_MASS3L, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SNOW_AGE,       & 
       SNOW_FLG3L,                                                      &
! Lestevens Sept 2012 - adding progs for CASACNP
       CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
       NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
       cable_lai,                                                       &

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC, &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC, &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC, &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC, TRACER_LBC, &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND, &
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND, &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND, &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND, &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND, &
       MURK_LBC_TEND, TRACER_LBC_TEND, &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB, &
      ! 2.2: Data variables stored in secondary space.
       PSTAR, &
      ! 2.3: Cloud fields
       CCB, CCT, CCLWP, &
      ! 2.4: Boundary layer fields
       ZH, & 
      ! Standard deviation of turbulent fluctuations of layer 1
       T1_SD, &
      ! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD, &
      ! Number of model levels in the  turbulently mixed layer
       NTML, &
      ! Top level for turb mixing in any decoupled Sc layer
       NTDSC, &
      ! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS, & 
      ! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT, &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN, &
      ! 2.5: Vegetation Ancillary fields
       CANOPY_WATER, SURF_CAP, SURF_RESIST, ROOT_DEPTH, INFILT, &
       VEG_FRAC, LAI, CANHT, Z0, SFA, MDSA, GS, &
      ! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2, &
       OROG_GRAD_X, OROG_GRAD_Y, &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY, &
      ! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS, &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, &
      ! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT, &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND, &
       DUST_MREL1, DUST_MREL2, DUST_MREL3, &
       DUST_MREL4, DUST_MREL5, DUST_MREL6, &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM, &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX, &
      ! Tracer Fluxes - kdcorbin, 05/10
       TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4, &
       TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8, &
       TRACER_FLUX9, TRACER_FLUX10, TRACER_FLUX11, TRACER_FLUX12, &
       TRACER_FLUX13, TRACER_FLUX14, TRACER_FLUX15, TRACER_FLUX16, &
       TRACER_FLUX17, TRACER_FLUX18, TRACER_FLUX19, TRACER_FLUX20, &

      ! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5, &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10, &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15, &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20, &
      !   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX, &
      !   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5, &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9, &
       LAI_PFT, CANHT_PFT, DISTURB_VEG, &
       SOIL_ALB, SOIL_CARB, &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4, &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC, &
       RSP_W_PFT_ACC, RSP_S_ACC, &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4, &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE, &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, &
       DOLR_FIELD, &
       LW_DOWN, SW_TILE_RTS, &
!! MODIFIED BY AT
!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!       TSLAB, TCLIM, HCLIM, CHEAT, OIFLX, UICE, VICE, &
!       SIG11NE, SIG11SE, SIG11SW, SIG11NW, &
!       SIG12NE, SIG12SE, SIG12SW, SIG12NW, &
!       SIG22NE, SIG22SE, SIG22SW, SIG22NW, &
!! END MODIFIED BY AT
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS, &
!   2.15: Fields carried forward from previous version.
!         May not be required
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho, &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP, &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND, &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET, &
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE, &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM, &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT, &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN, &
       C_SOLAR,C_BLUE,C_DOWN,C_LONGWAVE,C_TAUX,C_TAUY,C_WINDMIX, &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_BOTMELTN,C_TOPMELTN,C_LSRAIN, &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_PRESS, C_U10, C_V10, &
! UKCA oxidant fields
        OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA, & 
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR, &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK, &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4, &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC, &
       ARCLDLTA_DL, &
! Convective Cloud Fields
       LCBASE, CCW_RAD, &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM, &

! End arg_atm_fields.h
! ARGBND Control data calculated from NAMELIST-
     &  NBOUND_LOOKUP,                                                  &
      ! Headers from atmosphere boundary data sets
      ! Headers from ocean boundary data sets
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
! Constant arrays needed for Atmosphere-routing coupling
! ARGATCPL start
! Description: Gridline coordinates for interpolation and area-averaging
! between atmosphere and river-routing grids (Part of ARGAOCPL.hdk)
! Author: C.Bunton 28.02.03
!
! History:
! Version  Date    Comment
!  5.5  28/02/03  Original code. C.Bunton
     &  XPA, XUA, XVA, YPA, YUA, YVA,                                   &
! END ARGATCPL
     & G_P_FIELD,                                                       &
     & G_R_FIELD,                                                       &
     & OBS_FLAG,OBS,obs_flag_len,obs_len,                               &
     &      ngrgas,grgas_addr,                                          &
     &                     IAU_lookup,                                  &
                                       ! inout
     &                     D1_IAU_k4,                                   &
                                       ! inout
     &                     D1_TFilt,   &! inout
     &                     endstep)
      USE rad_switches_mod
      Use cv_cntl_mod, Only:                                            &
          lcv_phase_lim, lcv_3d_cca, lcv_3d_ccw, lcv_ccrad,             &
          lcv_pc2_diag_sh             
      Use cv_run_mod,  Only:                                            &
          l_rp, l_rp2, convection_option
      Use level_heights_mod
      Use trignometric_mod
      Use dyn_coriolis_mod
      Use dyn_var_res_mod
      Use diff_coeff_mod
      Use rad_mask_trop_mod
      Use rot_coeff_mod
      use cable_iday_mod
      IMPLICIT NONE
!
! Description: Perform a 1-timestep integration of the Atmosphere Model,
!   including assimilation, physics and dynamics processing.
!
! Method: Sequential execution of code sections below, with control and
!  logical switches determining whether each section is accessed.
!
!!ND ND code section and UM added, tidy and re-number required.
! Section 0.  Initialisation.
! Section 0.1  Filter near polar winds and theta
! Section 1.0  Call Atmospheric Physics1
!     STASH diagnostics for Physics1.
! Section 2. Calculate new time-level estimates to Theta, q, qcl, qcf
! Section 2.1 Calculate advective Momentum increments.
! Section 2.2 Calculate Semi-Lagrangian part of continuity equation.
! Section 2.3 Diagnostics at end of advection.
!     STASH diagnostics for Dynamics advection.
! Section 3.0  Call Atmospheric Physics2
!     STASH diagnostics for Physics2.
! Section 4.  Form and solve Helmholtz equation and update variables.
! Section 4.1 Form and solve Helmholtz equation, return all corrections
!     Set lateral boundaries.
! Section 5.0 Calculate rho at new time level, using flux form.
! Section 6.0 Calculate u, v, w, theta, q, qcl, qcf at new time level.
! Section 7.0 Mean all polar variables on a level to remove deviation.
! Section 8.0 Update pressure to value at new-time level.
!     VAR or AC assimilation.
! Section 9.0 Diagnostics at end of timestep
!     STASH diagnostics - dynamics based (section 15).
!     STASH diagnostics - physics  based (section 16).
!     STASH diagnostics - climate end of timestep (section 30).
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component: Control Atmos
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
! CMAXSIZE maximum sizes for dimensioning arrays
! of model constants whose sizes are configuration dependent. This
! allows constants to be read in from a NAMELIST file and maintain
! the flexibility of dynamic allocation for primary variables. The
! maximum sizes should agree with the maximum sizes implicit in the
! front-end User Interface.

!
!  Model            Modification history:
! version  Date
! 3.2  26/03/93  New COMDECK. Author R.Rawlins
! 3.4  06/08/94: Parameter MAX_NO_OF_SEGS used to dimension addresses
!                in macro-tasked calls to SWRAD, LWRAD & CONVECT.
!                Authors: A.Dickinson, D.Salmond, Reviewer: R.Barnes
! 3.5  22/05/95  Add MAX_N_INTF. D. Robinson
! 4.5  29/07/98  Increase MAX_N_INTF/MAX_N_INTF_A to 8. D. Robinson.
! 5.0  20/04/99  Changes for conversion to C-P C dynamics grid.
!                R. Rawlins
!  6.1   04/08/04  Add diffusion variable max_power     Terry Davies
! 6.2  25/12/05  Add max_updiff_levels/max_sponge_width   Terry Davies

      INTEGER,PARAMETER::max_model_levels = 100 ! Maximum no. of levels

      ! Max levels in boundary layer
      INTEGER,PARAMETER:: max_bl_levels = max_model_levels

      ! Max size of alpha_Cd
      INTEGER,PARAMETER :: max_number_alpha_cds = max_bl_levels

      ! Max no. of levels for pvort output
      INTEGER,PARAMETER :: MAX_REQ_THPV_LEVS = max_model_levels

      ! Max no. 1-2-1 rows in polar filter
      INTEGER,PARAMETER ::  max_121_rows =  8
      ! 0 is used for horizontal diffusion pointer

      ! Max no. of levels (from top) to apply upper level diffusion
      INTEGER,PARAMETER ::  max_updiff_levels = 10

      ! Max size of any sponge zones
      INTEGER,PARAMETER ::  max_sponge_width = 10

      ! Max size of look-up tables for searches
      INTEGER,PARAMETER ::  max_look = 2048

      ! Max no. of atmos interface areas
      INTEGER,PARAMETER :: MAX_N_INTF_A =  8

      ! Max no. of points in LBC      
      INTEGER,PARAMETER :: max_intf_lbcrow_length = 1000
      INTEGER,PARAMETER :: max_intf_lbcrows = 1000
        
      ! Max no. of atmos interface levels
      INTEGER,PARAMETER :: MAX_INTF_LEVELS = max_model_levels

      ! Maximum number of physics segments
      INTEGER,PARAMETER :: MAX_NO_OF_SEGS = 200
      ! MAX_N_INTF/MAX_N_INTF_A to be sorted out in later version
      ! Max no. of interface areas
      INTEGER, PARAMETER :: MAX_N_INTF =  8
! CMAXSIZE end
! CSUBMODL start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.5    07/04/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!  1. Internal model and submodel dump partition identifiers - fixed
!     for all experiments.
! CSMID start
!
! Description:
!    Hold parameters defining internal model identifiers and submodel
!    data partition (ie main D1 data array and consequent dump), both
!    short and long form.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.3    26/10/93   M. Carter. Part of an extensive mod that:
!                    1.Removes the limit on primary STASH item numbers.
!                    2.Removes the assumption that (section,item)
!                      defines the sub-model.
!                    3.Thus allows for user-prognostics.
!                    Add index to submodel home dump.
! 3.5    13/03/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
! 6.0    02/07/03   Add X_IM and X_SM for small exec.      E.Leung
!
! Declarations:
!
!   Hold parameters defining internal model identifiers and submodel
!   data partition (ie main D1 data array and consequent dump), both
!   short and long form
      ! Internal models
      INTEGER,PARAMETER:: A_IM      = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: ATMOS_IM  = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: O_IM      = 2 ! Ocean internal model
      INTEGER,PARAMETER:: OCEAN_IM  = 2 ! Ocean internalmodel
      INTEGER,PARAMETER:: S_IM      = 3 ! Slab internal model
      INTEGER,PARAMETER:: SLAB_IM   = 3 ! Slab internal model
      INTEGER,PARAMETER:: W_IM      = 4 ! Wave internal model
      INTEGER,PARAMETER:: WAVE_IM   = 4 ! Wave internal model
      INTEGER,PARAMETER:: I_IM      = 5 ! Sea=ice internal model
      INTEGER,PARAMETER:: SEAICE_IM = 5 ! Sea=ice internal model
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_IM      = 6 ! ND internal model
      INTEGER,PARAMETER:: NATMOS_IM = 6 ! ND internal model
      ! Small Executables
      INTEGER,PARAMETER:: X_IM      = 7 ! SX indicator

      ! Submodels
      INTEGER,PARAMETER:: A_SM      = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: ATMOS_SM  = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: O_SM      = 2 ! Ocean submodel
      INTEGER,PARAMETER:: OCEAN_SM  = 2 ! Ocean submodel
      INTEGER,PARAMETER:: W_SM      = 4 ! Wave submodel
      INTEGER,PARAMETER:: WAVE_SM   = 4 ! Wave submodel
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_SM      = 6 ! ND submodel
      INTEGER,PARAMETER:: NATMOS_SM = 6 ! ND submodel
      ! Small Executables
      INTEGER,PARAMETER:: X_SM      = 7 ! SX indicator

! CSMID end

!
!  2. Maximum internal model/submodel array sizes for this version.
!
! CSUBMAX start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.5    13/07/95   Original code. D.M. Goddard
! 4.0     3/11/95   Reduce max internal model, submodel from 10 to 4
!                   to save space in model. At 4.0 the max no of
!                   supported models is 3, 1 slot is reserved for
!                   expansion. Rick Rawlins.
!  4.1  21/02/96  Wave model introduced as 4th sub-model.  RTHBarnes
!
! Declarations:
!
!
!  1. Maximum internal model/submodel array sizes for this version.
!
      ! Max no. of internal models
      INTEGER,PARAMETER:: N_INTERNAL_MODEL_MAX=4

      ! Max no. of submodel dump partitions
      INTEGER,PARAMETER:: N_SUBMODEL_PARTITION_MAX=4

      ! Max value of internal model id
      INTEGER,PARAMETER:: INTERNAL_ID_MAX=N_INTERNAL_MODEL_MAX

      ! Max value of submodel dump id
      INTEGER,PARAMETER:: SUBMODEL_ID_MAX=N_SUBMODEL_PARTITION_MAX

! CSUBMAX end
!
!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.
      INTEGER :: N_INTERNAL_MODEL          ! No. of internal models
      INTEGER :: N_SUBMODEL_PARTITION      ! No. of submodel partitions

      ! Internal models
      INTEGER :: INTERNAL_MODEL_LIST(N_INTERNAL_MODEL_MAX)

      ! Submodel identifier for each internal model in list
      INTEGER :: SUBMODEL_FOR_IM    (N_INTERNAL_MODEL_MAX)

      ! Submodel number for each submodel id
      INTEGER :: SUBMODEL_FOR_SM(N_INTERNAL_MODEL_MAX)

      ! Namelist for information in 3.
      NAMELIST/NSUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,          &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM

      ! 4. Lists calculated in model from user interface supplied arrays
      ! experiment specific.

      ! No of internal models in each submodel partition indexed by sm
      !  identifier
      INTEGER :: N_INTERNAL_FOR_SM(SUBMODEL_ID_MAX)

      ! List of  submodel partition identifiers
      INTEGER :: SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION_MAX)

      ! Submodel partition identifier indexed by internal model identifie
      INTEGER :: SUBMODEL_PARTITION_INDEX(INTERNAL_ID_MAX)

      ! Sequence number of internal model indexed by internal model
      ! identifier: required to map from id to STASH internal model
      ! sequence
      INTEGER :: INTERNAL_MODEL_INDEX(INTERNAL_ID_MAX)


      ! Last internal model within a submodel partition if .TRUE.,
      ! indexed by internal model id.
      LOGICAL :: LAST_IM_IN_SM(INTERNAL_ID_MAX)

      ! Common block for information in 3. and 4.
      COMMON/SUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,             &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM,SUBMODEL_FOR_SM,            &
     &  N_INTERNAL_FOR_SM,SUBMODEL_PARTITION_LIST,                      &
     &  SUBMODEL_PARTITION_INDEX,                                       &
     &  INTERNAL_MODEL_INDEX,                                           &
     &  LAST_IM_IN_SM

!
!  5. Time information specifying coupling frequencies between internal
!     models and submodels, and multipliers, indexed by sequence of
!     internal models and submodels (ie left to right along node tree).
!     {Not required at this release}.
!
! Namelists for information in 5. {Not required at this release}
!
!
!  6. Lists of coupling nodes defining coupling frequencies between
!     internal models and between submodel partitions. (Not defined
!     yet at this release).
!CALL CNODE
!
!  7. Variables dealing with general coupling switches at the control
!     level. {These will require revision at the next release when
!     coupling between internal models is dealt with more generally.
!     Logicals below are set in routine SETGRCTL.}

      ! new internal model next group of timesteps if .true.
      LOGICAL :: new_im

      ! new submodel dump  next group of timesteps if .true.
      LOGICAL :: new_sm

      COMMON/CSUBMGRP/new_im,new_sm

      INTEGER SUBMODEL_IDENT
      COMMON/SUBMODID/SUBMODEL_IDENT
! CSUBMODL end
! ------------------------ Comdeck PARVARS -------------------------
! Parameters and common blocks required by the mpp-UM
! Added new comdeck AMAXSIZE required for new arrays in PARCOMM
! Add non-mpp option
!                                                      P.Burton
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end
!====================== COMDECK AMAXSIZE ========================
! Description
!   This comdeck provides parameters giving the maximum likely sizes
!   of key UM resolution variables, useful for sizing static arrays.
!
!   History:
!   Model    Date     Modification history
!  version
!   4.2      18/11/96 New comdeck created.  P.Burton
!   4.3      24/01/97 Define MaxFieldSize to be a quarter of the
!                     SHMEM common block size.         P.Burton
!   4.4      3/7/97   Add MaxFieldSizeMes. Deborah Salmond
!   4.5     12/01/98  Added new variables, and changed sizes to
!                     correspond to global hi-res forecast - current
!                     largest configuration.                P.Burton
!                     Changed MAX_SHMEM_COMMON_SIZE to 3000000
!                     required for operational data assimilation.
!                                                           P.Burton
!   5.0     29/04/99  Changed variable names:
!                       P_ROWS_MAX -> ROWS_MAX
!                       P_LEVELS_MAX -> MODEL_LEVELS_MAX
!                       Q_LEVELS_MAX -> WET_LEVELS_MAX
!                       MaxHaloSize -> MaxHaloArea
!                     Removed variable:
!                       HALO_MAX (use PARPARM Max_Halo_Size instead)
!    5.0   29/04/99  Remove mpp #define
!    5.3   05/12/01  Remove MaxFieldSize, MaxFieldSizeMes and
!                    Max3DFieldSize.  S.Cusack
!    5.5   22/01/03  Increase ROW_LENGTH_MAX and HORIZ_DIM_MAX
!                    from 432 to 548. D Robinson.
!    6.1   31/08/04  Allow up to 100 levels.  R.Barnes
!    6.2   13/02/06  Increase max values of row_length and
!                    rows to cope with FOAM high res, as well
!                    as Global N320 and NAE.  M Martin.
!    6.2   24/11/05  Use max function for horiz_dim_max. R Barnes
!    6.2     11/01/06 Remove max_shmem_common_size here and
!                     in rdobs2.   Camilla Mathison/R Barnes
!

! Maximum sector size for I/O
      INTEGER,PARAMETER:: IO_SECTOR_SIZE_MAX=4096
      INTEGER,PARAMETER:: ROW_LENGTH_MAX   = 840 ! Maximum row length
      INTEGER,PARAMETER:: ROWS_MAX         = 600 ! Max no of rows

      ! MAX(ROW_LENGTH_MAX,ROWS_MAX)
      INTEGER,PARAMETER:: HORIZ_DIM_MAX=MAX(ROW_LENGTH_MAX,ROWS_MAX)

      INTEGER,PARAMETER:: MODEL_LEVELS_MAX = 100 ! Max no of total levels
      INTEGER,PARAMETER:: WET_LEVELS_MAX   = 100 ! Max no of wet levels
      INTEGER, PARAMETER :: Max2DFieldSize = ROW_LENGTH_MAX*ROWS_MAX +  &
     &  IO_SECTOR_SIZE_MAX
      INTEGER, PARAMETER :: MaxHaloArea    = HORIZ_DIM_MAX*Max_Halo_Size
!========================== COMDECK PARCOMM ====================
!
! *** NOTE : This comdeck requires comdeck PARPARM to be *CALLed
!            first.
!
!   Description:
!
!   This COMDECK contains COMMON blocks for the mpp-UM
!
!
!   Two COMMON blocks are defined:
!     i)  UM_PARVAR holds information required by the
!         Parallel Unified Model itself
!     ii) MP_PARVAR holds information required by the interface to
!         the Message Passing Software used by the PUM
!
!   Key concepts used in the inline documentation are:
!     o global data - the entire data domain processed by the UM
!     o LOCAL data - the fragment of the global data which is
!       stored by this particular process
!     o PERSONAL data - the fragment of the LOCAL data which is
!       updated by this particular process
!     o HALO data - a halo around the PERSONAL data which forms
!       the LOCAL data
!
!     Acronyms used:
!     LPG - Logical Process Grid, this is the grid of logical
!           processors; each logical processor handles one of the
!           decomposed parts of the global data. It does not
!           necessarily represent a physical grid of processors.
!
!   History:
!
!   4.1      27/1/96  New comdeck based on second section of
!                     old PARVARS.   P.Burton
!   4.2     19/08/96  Removed some unused variables, and added
!                     current_decomp_type variable to allow use
!                     of flexible decompositions.
!                     Added nproc_max to indicate the max. number
!                     of processors used for mpp-UM
!                                                      P.Burton
!   5.0     12/04/99  - Added halosize array to common block
!                     - Added halo_i and halo_j to common block
!                     - Added fld_type dimension to glsize
!                     - Added halo type dimension to lasize
!                     - Added fld_type dimension to lasize
!                     - Replaced blsizep/blsizeu by blsize with
!                       extra fld_type dimension
!                     - Replace attop etc. with at_extremity
!                     - Added g_pe_index to common block
!                                                      P.Burton
!   5.1     22/05/00  Removed DATA statement and put in BLKDATA
!                                                      P.Burton
!   5.1     26/01/00  - Renamed g_pe_index -> g_pe_index_EW
!                     - Added g_pe_index_NS
!                                                     P.Burton
!   5.2     02/08/00  Added g_at_extremity        P.Burton
!   5.3     14/09/01  Added sb_model_domain   P.Burton
!   5.5     06/08/00  Modification for parallelisation of WAM.
!                     Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!   5.5     30/01/03  Generalised datastart   P.Selwood.
!   5.5     07/02/03  SX now uses PARCOMM instead of SXCOMM    E.Leung
!   6.0     18/09/03  F90-fy continuation lines.               P.Dando
!   6.2     23/11/05  Removed all references to the wavemodel.
!                     T.Edwards
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for the Parallel Unified Model
! =======================================================

      INTEGER :: first_comp_pe       ! top left pe in LPG
      INTEGER :: last_comp_pe        ! bottom right pe in LPG
      INTEGER :: current_decomp_type ! current decomposition type
      INTEGER :: Offx                ! standard halo size in East-West
      INTEGER :: Offy                ! standard halo size in North-South
      INTEGER :: halo_i              ! extended halo size in East-West
      INTEGER :: halo_j              ! extended halo size in North-South
      INTEGER :: halosize(Ndim_max,NHalo_max) ! available halo sizes
      INTEGER :: glsize(Ndim_max,Nfld_max) ! global data size
      INTEGER :: lasize(Ndim_max,Nfld_max,NHalo_max) ! local data size
      INTEGER :: blsize(Ndim_max,Nfld_max) ! personal data size

      ! position of personal data in global data (in terms of standard
      ! Fortran array notation
      INTEGER :: datastart(Ndim_max)

      ! Generalised version of datastart for *all* fieldtypes
      INTEGER :: datastart_f(Ndim_max,Nfld_max)

      INTEGER :: gridsize(Ndim_max)  ! size of the LPG in each dimension

      ! position of this process in the LPG 0,1,2,...,nproc_x-1 etc.
      INTEGER :: gridpos(Ndim_max)

      INTEGER :: sb_model_domain


      ! logicals indicating if a processor is at the edge of the LPG
      LOGICAL :: at_extremity(4)

      COMMON /UM_PARVAR/                                                &
     &  first_comp_pe, last_comp_pe, current_decomp_type, Offx, Offy,   &
     &  halo_i, halo_j, halosize, glsize, lasize, blsize, datastart,    &
     &  datastart_f, gridsize, gridpos                                  &
     &,                 at_extremity,sb_model_domain

      ! Common block for the Message Passing Software

      ! type of boundary (cyclic or static) in each direction
      INTEGER :: bound(Ndim_max)

      ! global copy of local data size
      INTEGER :: g_lasize(Ndim_max,Nfld_max,NHalo_max,0:maxproc)

      ! global copy of personal data area
      INTEGER :: g_blsize(Ndim_max,Nfld_max,0:maxproc)

      ! global copy of datastart
      INTEGER :: g_datastart(Ndim_max,0:maxproc)

      ! global copy of datastart_f
      INTEGER :: g_datastart_f(Ndim_max,Nfld_max,0:maxproc)

      INTEGER :: g_gridpos(Ndim_max,0:maxproc) ! global copy of gridpos

      ! Which processor column a given point is in: 0 -> nproc_x-1
      INTEGER :: g_pe_index_EW(1-Max_Halo_Size:                         &
     &  ROW_LENGTH_MAX+Max_Halo_Size)

      ! Which processor row a given point is in: 0 -> nproc_y-1
      INTEGER :: g_pe_index_NS(1-Max_Halo_Size:ROWS_MAX+Max_Halo_Size)

      INTEGER :: nproc      ! number of processors in current decomp
      INTEGER :: mype      ! number of this processor (starting from 0)
      INTEGER :: nproc_max  ! maximum number of processors
      INTEGER :: nproc_x    ! number of processors in x-direction
      INTEGER :: nproc_y    ! number of processors in y-direction

      ! array with the tids of the four neighbours in the horizontal
      ! plane
      INTEGER :: neighbour(4)

      INTEGER :: gc_proc_row_group  ! GID for procs along a proc row
      INTEGER :: gc_proc_col_group  ! GID for procs along a proc col
      INTEGER :: gc_all_proc_group  ! GID for all procs

      ! at_extremity for each processor
      LOGICAL :: g_at_extremity(4,0:maxproc)

      COMMON /MP_PARVAR/                                                &
     &  bound,                                                          &
     &  g_lasize,g_blsize,g_datastart, g_datastart_f, g_gridpos,        &
     &  g_pe_index_EW,g_pe_index_NS,                                    &
     &  nproc_max,nproc_x,nproc_y,                                      &
     &  neighbour,gc_proc_row_group,                                    &
     &  gc_proc_col_group, gc_all_proc_group                            &
     &  ,nproc,mype                                                     &
     &, g_at_extremity

! PARCOMM end

! TYPSIZE start
!   Description:
!     This file contains sizes needed for dynamic allocation of
!   main data arrays within the model. Sizes read in from the user
!   interface via NAMELISTs are passed by /COMMON/. Other control
!   sizes that are fundamental in the definition of data structures
!   are assigned by PARAMETER statements.
!
      ! All sizes
      ! Not dependent on sub-model
      ! DATA IN NAMLST#x MEMBER OF THE JOB LIBRARY
      ! atmos START
      ! Main sizes of fields for each submodel
      ! Grid-related sizes for ATMOSPHERE submodel.
      INTEGER:: ROW_LENGTH           ! IN: No of points per local row
      INTEGER:: global_ROW_LENGTH    ! IN: Points per global row
      INTEGER:: ROWS                 ! IN: No of local (theta) rows
      INTEGER:: global_ROWS          ! IN: No of global (theta) rows
      INTEGER:: MODEL_LEVELS         ! IN: No of model levels
      INTEGER:: LAND_FIELD           ! IN: No of land points in field
      INTEGER:: NTILES               ! IN: No of land surface tiles

      ! Physics-related sizes for ATMOSPHERE submodel
      INTEGER:: WET_LEVELS          ! IN: No of moist-levels
      INTEGER:: CLOUD_LEVELS        ! IN: No of cloud-levels
      INTEGER:: ST_LEVELS           ! IN: No of soil temperature levels
      INTEGER:: SM_LEVELS           ! IN: No of soil moisture levels
      INTEGER:: BL_LEVELS           ! IN: No of boundary-layer-levels
      INTEGER :: OZONE_LEVELS       ! IN: No of ozone-levels
      INTEGER :: TPPS_OZONE_LEVELS  ! IN: No of tropopause-ozone-levels
      INTEGER :: RIVER_ROWS         ! IN: No of rows for river routing
      INTEGER :: RIVER_ROW_LENGTH   ! IN: Row length for river routing
      ! Dynamics-related sizes for ATMOSPHERE submodel

      INTEGER:: TR_LEVELS            ! IN: No of tracer-levels
      INTEGER:: TR_VARS              ! IN: No of passive tracers
      INTEGER:: TR_UKCA              ! IN: No of UKCA tracers

      ! Dynamics output diagnostic-related sizes for ATMOSPHERE submodel
      INTEGER:: THETA_PV_P_LEVS   ! IN: No of levels requested for pvort

      ! For Small executables
      INTEGER:: TOT_LEVELS
      ! Assimilation-related sizes for ATMOSPHERE submodel
      INTEGER :: N_AOBS           ! IN: No. of atmos observation types

      ! Grid related sizes for data structure
      ! Data structure sizes for ATMOSPHERE submodel
      INTEGER:: A_PROG_LOOKUP     ! IN: No of prognostic fields
      INTEGER:: A_PROG_LEN        ! IN: Total length of prog fields
      INTEGER:: A_LEN_INTHD       ! IN: Length of INTEGER header
      INTEGER:: A_LEN_REALHD      ! IN: Length of REAL header
      INTEGER:: A_LEN2_LEVDEPC    ! IN: No of LEVEL-dependent arrays
      INTEGER:: A_LEN2_ROWDEPC    ! IN: No of ROW-dependent arrays
      INTEGER:: A_LEN2_COLDEPC    ! IN: No of COLUMN-dependent arrays
      INTEGER:: A_LEN2_FLDDEPC    ! IN: No of FIELD arrays
      INTEGER:: A_LEN_EXTCNST     ! IN: No of EXTRA scalar constants
      INTEGER:: A_LEN_CFI1        ! IN: Length of compressed fld index 1
      INTEGER:: A_LEN_CFI2        ! IN: Length of compressed fld index 2
      INTEGER:: A_LEN_CFI3        ! IN: Length of compressed fld index 3
      ! atmos end

      ! OCEAN start
! TYPOCPAR start
!  History:
!  Version   Date     Comment
!  -------   ----     -------
!    4.4   15.06.97   Add free surface scalar R.Lenton
!     5.1   07.01.00   Invert_ocean to false with New Dynamics. JC Thil.
!     5.4   29.08.02   Add N_STRAIT and N_STRAIT_CLM. D. Storkey
!     5.5   03.01.03   Remove typocbas. R. Hill
!
      ! Grid related sizes for OCEAN model
      INTEGER::LSEG!IN:Max no of sets of start/end indices for vorticity
      INTEGER:: NISLE            ! IN: No of islands
      INTEGER:: ISEGM            ! IN: Max no of island segments per box
      INTEGER :: N_STRAIT       ! IN: No of pairs of Strait exchange pts
      INTEGER :: N_STRAIT_CLM   ! IN: No of Strait pts set by climate
                                !     values
      INTEGER:: O_LEN_COMPRESSED ! IN: No of ocean points in 3D field
      INTEGER:: LSEGC            ! IN: No of island basins for mead calc
      ! No of start/end indicies for the free surface solution
      INTEGER:: LSEGFS ! IN

      ! Fourier filtering for OCEAN submodel
      INTEGER :: LSEGF    ! IN: max. no of sets of indices for filtering
      INTEGER :: JFRST    ! IN: first J row of T to be filtered

      ! filtering is done on T with a low pass cut off to make the zonal
      ! dimension of the box filtered effectively the same as that of
      ! the boxes on row JFT0
      INTEGER :: JFT0     ! IN:

      INTEGER :: JFT1     ! IN: last J row of T in SH to be filtered
      INTEGER :: JFT2     ! IN: first J row of T in NH to be filtered
      INTEGER :: JFU0     ! IN: same function as JFT0 but for U,V
      INTEGER :: JFU1     ! IN: last J row of U,V in SH to be filtered
      INTEGER :: JFU2     ! IN: first J row of U,V in NH to be filtered

      ! Variables derived from those above
      INTEGER :: IMU     ! IN: total number of U,V grid boxes zonally
      INTEGER :: IMTP1   ! IN: IMT+1
      INTEGER :: IMTM1   ! IN: IMT-1
      INTEGER :: IMTM2   ! IN: IMT-2
      INTEGER :: IMUM1   ! IN: IMU-1
      INTEGER :: IMUM2   ! IN: IMU-2
      INTEGER :: JMTP1   ! IN: JMT+1
      INTEGER :: JMTM1   ! IN: JMT-1
      INTEGER :: JMTM2   ! IN: JMT-2
      INTEGER :: JSCAN   ! IN: JMTM2+1
      INTEGER :: KMP1    ! IN: KM+1
      INTEGER :: KMP2    ! IN: KM+2
      INTEGER :: KMM1    ! IN: KM-1
      INTEGER :: NSLAB   ! IN: no of words in one slab
      INTEGER :: JSKPT   ! IN: no of rows of T and U,V not filtered in
      INTEGER :: JSKPU   ! IN: low and mid latitudes + 1
      INTEGER :: NJTBFT  ! IN: no of J rows to be filtered on T
      INTEGER :: NJTBFU  ! IN: no of J rows to be filtered on U,V
      INTEGER :: IMTKM   ! IN: IMT*KM
      INTEGER :: NTMIN2  ! IN: maximum of NT or 2
      INTEGER :: IMTD2   ! IN: IMT/2
      INTEGER :: LQMSUM  ! IN: IMTD2*(IMT-IMTD2)
      INTEGER :: LHSUM   ! IN: IMT*IMTP1/2
      INTEGER :: IMTX8   ! IN: IMT*8
      INTEGER :: IMTIMT  ! IN: IMT*IMT
      INTEGER :: IMROT   ! X dimension for Coriolis array
      INTEGER :: JMROT   ! Y dimension for Coriolis array
      INTEGER :: IMBC    ! No of columns in boundary field array
      INTEGER :: JMBC    ! No of rows in boundary field array
      INTEGER :: KMBC    ! No of levels in boundary field array
      INTEGER :: NTBC    ! No of tracers in boundary field array
      INTEGER :: JMMD    ! No of rows for mead diagnostic basin indices
      INTEGER :: LDIV    ! No of divisions mead basin indices

      ! Grid-related switches for OCEAN submodel
      LOGICAL :: CYCLIC_OCEAN        ! IN: TRUE if CYCLIC E-W boundary
      LOGICAL :: GLOBAL_OCEAN        ! IN: TRUE if global domain

      ! TRUE if ocean grid NS-inverted cf atmos
      LOGICAL, PARAMETER :: INVERT_OCEAN=.FALSE.
      ! User interface limit for tracers
      ! Max no. tracers in STASHMASTER
      INTEGER, PARAMETER :: O_MAX_TRACERS=20
! COMOCPAR start
      COMMON /COMOCPAR/ GLOBAL_OCEAN, CYCLIC_OCEAN,                     &
     &  LSEG,NISLE,ISEGM,N_STRAIT,N_STRAIT_CLM,O_LEN_COMPRESSED,LSEGC,  &
     &  LSEGFS,LSEGF,JFRST,JFT0,JFT1,JFT2,JFU0,JFU1,JFU2,IMU,IMTP1,     &
     &  IMTM1,IMTM2,IMUM1,IMUM2,JMTP1,JMTM1,JMTM2,JSCAN,KMP1,KMP2,KMM1, &
     &  NSLAB,JSKPT,JSKPU,NJTBFT,NJTBFU,IMTKM,NTMIN2,                   &
     &  IMTD2,LQMSUM,LHSUM,IMTX8,IMTIMT,                                &
     &  IMROT,JMROT,IMBC,JMBC,KMBC,NTBC,JMMD,LDIV
! COMOCPAR end
! TYPOCPAR end
! TYPOCBAS Physics-related sizes for OCEAN submodel
      INTEGER ::  NT ! IN: No of ocean tracers (inc T,S)
      ! Grid related sizes for OCEAN model
      INTEGER :: IMT  ! IN: No of points per row (incl wrap)
      INTEGER :: JMT  ! IN: No of tracer rows
      INTEGER :: KM   ! IN: No of tracer levels
      INTEGER :: NT_UI     ! Copy of NT
      INTEGER :: IMT_UI    ! Copy of IMT
      INTEGER :: JMT_UI    ! Copy of JMT
      INTEGER :: KM_UI     ! Copy of KM
      INTEGER :: NICE      ! IN: No. of sea ice thickness categories
! COMOCBAS start
      COMMON /COMOCBAS/ NT_UI, IMT_UI, JMT_UI, KM_UI                    &
     &        ,NT, IMT, JMT, KM                                         &
     &                 ,NICE
! COMOCBAS end
! TYPOCBAS end
! TYPOASZ sizes for dynamic allocation of ocean assim.
! 5.2 11/08/00  JO_NMAX_OBS_ICE introduced. JO_MAX_OBS_VAL and
!               JO_NMAX_OBS_ICE put into COMOCASZ. M J Bell
      INTEGER :: JO_MAX_OBS_VAL !max number of values in OBS array

      ! max no of flds reqd at once in sea ice assim
      INTEGER :: JO_NMAX_OBS_ICE

      ! length of climate/covariances array
      INTEGER, PARAMETER:: JO_LEN_COV = 1

      ! max number of columns in climate grid
      INTEGER,PARAMETER:: JO_MAX_COLS_C = 1

      ! max number of rows in climate grid
      INTEGER,PARAMETER:: JO_MAX_ROWS_C = 1

      ! max number of levels in climate grid
      INTEGER,PARAMETER:: JO_MAX_LEVS_C = 1

! COMOCASZ start
      COMMON /COMOCASZ/ JO_MAX_OBS_VAL, JO_NMAX_OBS_ICE
! COMOCASZ end
! TYPOASZ end
      ! OCEAN end

      !  WAVE sub-model start
      ! WAVE end

      ! Grid related sizes for COUPLING between atmos and OCEAN
      ! submodels [For mpp, sizes are 'global' values over all
      ! PEs.Also needed for river routing]
      INTEGER:: AOCPL_IMT                ! Ocean rowlength
      INTEGER:: AOCPL_JMT                ! Ocean no. of rows
      INTEGER:: AOCPL_ROW_LENGTH         ! Atmos rowlength
      INTEGER:: AOCPL_P_ROWS             ! Atmos no. of p rows

      COMMON/SIZE_AOCPL/                                                &
     &  AOCPL_IMT, AOCPL_JMT, AOCPL_ROW_LENGTH, AOCPL_P_ROWS

      ! Data structure sizes for ATMOSPHERE ANCILLARY file control
      ! routines
      INTEGER :: NANCIL_LOOKUPSA  ! IN: Max no of fields to be read

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER::N_INTF_A          ! No of atmosphere interface areas
      INTEGER::MAX_INTF_MODEL_LEVELS ! Max no of model levs in all areas
      INTEGER::MAX_LBCROW_LENGTH ! Max no of lbc row length in all areas
      INTEGER::MAX_LBCROWS ! Max no of lbc rows in all areas
      INTEGER::TOT_LEN_INTFA_P   ! Total length of interface p grids.
      INTEGER::TOT_LEN_INTFA_U    ! Total length of interface u grids.
      INTEGER::U_FIELD_INTFA      ! Length of Model U field (= U_FIELD)

      !  Data structure sizes for ATMOSPHERE BOUNDARY file control
      ! routines
! PARAMETERs defining the RIM_TYPE characteristics

!  History:
!  Date      Vn     Modification
!  31/10/01  5.3    Reset Nrima_max to 1. Change rima_type_orog
!                   from 2 to 1. D. Robinson

      INTEGER, PARAMETER :: Nrima_max = 1

      INTEGER, PARAMETER :: rima_type_norm=1  ! Normal field
      INTEGER, PARAMETER :: rima_type_orog=1  ! Orography field

! At 5.3 rima_type_orog=2 => rima_type_orog=1. This means that
! Orography LBCs has the same rim type as all prognostic LBCs. The
! value was changed rather than remove all occurences from the code to
! retain the functionality if it ever needs to be restored and to
! simplify the number of changes required.
      INTEGER :: RIMWIDTHA(Nrima_max)
      INTEGER :: NRIM_TIMESA      ! IN: Max no of timelevels in rim flds

      ! Data structure sizes for atmos & OCEAN BOUNDARY file control
      !routines

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)

      INTEGER :: PP_LEN_INTHD   ! IN: Length of PP file integer header
      INTEGER :: PP_LEN_REALHD  ! IN: Length of PP file real    header

      ! Other sizes passed from namelist into common blocks
      COMMON/NLSIZES/                                                   &
     &  ROW_LENGTH,global_ROW_LENGTH,ROWS,global_ROWS,                  &
     &  LAND_FIELD,MODEL_LEVELS,WET_LEVELS,                             &
     &  NTILES,                                                         &
     &  CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,           &
     &  OZONE_LEVELS,TPPS_OZONE_LEVELS,TR_VARS,TR_UKCA,                 &
     &  RIVER_ROWS, RIVER_ROW_LENGTH,                                   &

     &  THETA_PV_P_LEVS, N_AOBS,                                        &

     &  A_PROG_LOOKUP,A_PROG_LEN,                                       &
     &  A_LEN_INTHD,A_LEN_REALHD,                                       &
     &  A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,                   &
     &  A_LEN2_FLDDEPC,A_LEN_EXTCNST,                                   &
     &  A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,                               &

     &  NANCIL_LOOKUPSA,                                                &

     &  N_INTF_A, MAX_INTF_MODEL_LEVELS, MAX_LBCROW_LENGTH,             &
     &  MAX_LBCROWS, TOT_LEN_INTFA_P,                                   &
     &  TOT_LEN_INTFA_U, U_FIELD_INTFA,                                 &

     &  RIMWIDTHA, NRIM_TIMESA,                                         &

     &  PP_LEN_INTHD,PP_LEN_REALHD

      !-----------------------------------------------------------------
      ! data in STASHC#x member of the job library

      ! Data structure sizes for ATMOSPHERE submodel (config dependent)
      INTEGER:: A_LEN2_LOOKUP   ! IN: Total no of fields (incl diags)
      INTEGER:: A_LEN_DATA      ! IN: Total no of words of data
      INTEGER:: A_LEN_D1        ! IN: Total no of words in atmos D1
      ! Data structure sizes for SLAB submodel (config dependent)
      INTEGER:: S_LEN2_LOOKUP   !IN: Tot no of fields (incl diags)
      INTEGER:: S_LEN_DATA      !IN: Tot no of words of data
      ! Data structure sizes for OCEAN submodel (config dependent)
      INTEGER:: O_LEN2_LOOKUP     ! IN: Total no of fields (incl diags)
      INTEGER:: O_LEN_DATA        ! IN: Total no of words of data
      INTEGER:: O_LEN_DUALDATA    ! IN: Words of data at 2 time levels
      INTEGER:: O_LEN_D1          ! IN: Total no of words in ocean D1
      ! Data structure sizes for WAVE submodel (config dependent)
      INTEGER:: W_LEN2_LOOKUP     ! IN: Total no of fields (incl diags)
      INTEGER:: W_LEN_DATA        ! IN: Total no of words of data
      INTEGER:: W_LEN_D1          ! IN: Total no of words in atmos D1

      ! Size of main data array for this configuration

      INTEGER:: LEN_TOT             ! IN: Length of D1 array
      INTEGER:: N_OBJ_D1_MAX         ! IN: No of objects in D1 array

      COMMON/STSIZES/                                                   &
     &  S_LEN2_LOOKUP,S_LEN_DATA,                                       &
     &  A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                              &
     &  O_LEN2_LOOKUP,O_LEN_DATA,O_LEN_DUALDATA,O_LEN_D1,               &
     &  W_LEN2_LOOKUP,W_LEN_DATA,W_LEN_D1,                              &
     &  LEN_TOT,N_OBJ_D1_MAX
      ! global (ie. dump version) of *_LEN_DATA
      INTEGER:: global_A_LEN_DATA
      INTEGER:: global_O_LEN_DATA
      INTEGER :: global_W_LEN_DATA ! global (ie. dump version) of
                                   !                        W_LEN_DATA

      COMMON /MPP_STSIZES_extra/ global_A_LEN_DATA,global_O_LEN_DATA    &
     &      ,global_W_LEN_DATA
      ! Sizes of Stash Auxillary Arrays and associated index arrays
      ! Initialised in UMINDEX and UMINDEX_A/O/W
      INTEGER:: LEN_A_IXSTS
      INTEGER:: LEN_A_SPSTS
      INTEGER:: LEN_O_IXSTS
      INTEGER:: LEN_O_SPSTS
      INTEGER:: LEN_W_IXSTS
      INTEGER:: LEN_W_SPSTS

      COMMON /DSIZE_STS/                                                &
     &  LEN_A_IXSTS, LEN_A_SPSTS,LEN_O_IXSTS, LEN_O_SPSTS,              &
     &  LEN_W_IXSTS, LEN_W_SPSTS
!     From 4.5, the number of land points is computed for each
!     PE before the addressing section. All prognostics on land
!     points in the D1 space are now dimensioned by the local
!     no of land points rather than the global no of land points.

      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field
      ! ----------------------------------------------------------------
      ! extra variables not passed through user interface

      ! fundamental data sizes :
      ! Fundamental parameter  sizes of data structure
      ! Sizes applicable to all configurations (HISTORY FILE)

      ! Length of history file in dump
      INTEGER, PARAMETER :: LEN_DUMPHIST = 0

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      ! Length of dump fixed header
      INTEGER, PARAMETER :: LEN_FIXHD = 256

      ! Size of a single LOOKUP header
      INTEGER, PARAMETER :: LEN1_LOOKUP  = 64
      INTEGER, PARAMETER :: MPP_LEN1_LOOKUP= 2

      ! Size of compressed LBC LOOKUP (only used internally and
      ! contains just the items which change between each set of LBCs
      INTEGER, PARAMETER :: LEN1_LBC_COMP_LOOKUP = 8

      ! Sizes applicable to all configurations (STASH)
      ! Moved to <typstsz/typstsz.h>

      INTEGER:: INTF_LEN2_LEVDEPC !1st dim of interface out lev dep cons
      INTEGER:: INTF_LEN2_ROWDEPC !2nd dim of interface out Row dep cons
      INTEGER:: INTF_LEN2_COLDEPC !2nd dim of interface out Col dep cons
      
      COMMON /DSIZE/                                                    &
     &  INTF_LEN2_LEVDEPC,INTF_LEN2_ROWDEPC,INTF_LEN2_COLDEPC
      ! sub-model atmosphere   :
      ! Data structure sizes derived from grid size
      INTEGER:: A_LEN1_LEVDEPC ! IN: 1st dim of level  dep const
      INTEGER:: A_LEN1_ROWDEPC ! IN: 1st dim of row    dep const
      INTEGER:: A_LEN1_COLDEPC ! IN: 1st dim of column dep const
      INTEGER:: A_LEN1_FLDDEPC ! IN: 1st dim of field  dep const

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER :: INTF_LOOKUPSA        ! No of interface lookups.
      INTEGER :: OLD_INTF_LOOKUPSA    ! No of interface lookups
                                      ! for old (4.5) LBCs.
      COMMON /DSIZE_A/                                                  &
     &  A_LEN1_LEVDEPC,A_LEN1_FLDDEPC,A_LEN1_ROWDEPC,A_LEN1_COLDEPC,    &
     &  INTF_LOOKUPSA,OLD_INTF_LOOKUPSA

      ! sub-model atmosphere   : derived sizes
      ! derived from model grid/levels. Arakawa B-grid

                                  ! Size of fields on THETA grid:
      INTEGER :: THETA_FIELD_SIZE     ! IN: with no halos
      INTEGER :: THETA_OFF_SIZE       ! IN: with simple halos
      INTEGER :: THETA_HALO_SIZE      ! IN: with extended halos

                                  ! Size of fields on U grid:
      INTEGER :: U_FIELD_SIZE         ! IN: with no halos
      INTEGER :: U_OFF_SIZE           ! IN: with simple halos
      INTEGER :: U_HALO_SIZE          ! IN: with extended halos

                                  ! Size of fields on V grid
      INTEGER :: V_FIELD_SIZE         ! IN: with no halos
      INTEGER :: V_OFF_SIZE           ! IN: with simple halos
      INTEGER :: V_HALO_SIZE          ! IN: with extended halos

      INTEGER :: N_ROWS               ! IN: No of V-rows
      INTEGER :: N_CCA_LEV            ! IN: No of CCA Levels
      COMMON/DRSIZE_A/                                                  &
     &  N_ROWS,N_CCA_LEV,THETA_FIELD_SIZE,U_FIELD_SIZE,V_FIELD_SIZE,    &
     &  THETA_OFF_SIZE,U_OFF_SIZE,V_OFF_SIZE,                           &
     &  THETA_HALO_SIZE,U_HALO_SIZE,V_HALO_SIZE
      ! boundary updating      : derived values
      ! Variables describing the Atmosphere Lateral Boundary Conditions
      ! Local (per processor) information

      ! Size of atmos LBC for given field type, halo type and rimwidth
      ! type
      INTEGER:: LENRIMA(Nfld_max,NHalo_max,Nrima_max)

      ! Size of given side (PNorth,PEast,PSouth and PWest), field type,
                        ! halo type and rimwidth type
      INTEGER:: LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max)

                        ! Start of a given side within the LBC
      INTEGER:: LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max)

      ! Start of a given side within the LBC on a given processor
      INTEGER:: g_LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max,0:Maxproc-1)

      ! and global (within the file) information

      ! Size of atmos LBC on disk for given field type, halo type and
      ! rimwidth type
      INTEGER:: global_LENRIMA(Nfld_max,NHalo_max,Nrima_max)

                        ! Size of given side, field type and halo type
      INTEGER:: global_LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max)

                        ! Start of a given side within the LBC
      INTEGER:: global_LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max)

      ! Variables describing the Ocean Lateral Boundary Conditions
      INTEGER:: LENRIMO                ! Size of ocean LBC (theta)
      INTEGER:: LENRIMO_U              ! Size of ocean LBC (velocity)

      ! Variables that may be needed for vn5.2 but have not yet been
      ! dealt with at vn5.1
      INTEGER:: RIMFLDSA
      INTEGER:: RIMFLDSO
      INTEGER:: global_LENRIMDATA_A
      INTEGER:: LENRIMDATA_A
      INTEGER:: LENRIMDATA_O
      INTEGER:: BOUNDFLDS
      INTEGER:: RIM_LOOKUPSA
      INTEGER:: RIM_LOOKUPSO
      INTEGER:: BOUND_LOOKUPSA
      INTEGER:: BOUND_LOOKUPSO
      INTEGER:: BOUND_LOOKUPSW
      INTEGER :: RIM_LOOKUPSW
      INTEGER :: LENRIMDATA_W
      INTEGER :: global_LENRIMDATA_W
      COMMON/DRSIZ_BO/                                                  &
      ! Atmosphere variables
     &  LENRIMA, LBC_SIZEA, LBC_STARTA, g_LBC_STARTA,                   &
     &  global_LENRIMA,global_LBC_SIZEA,global_LBC_STARTA,              &
      ! Wave model variables
     & RIM_LOOKUPSW, LENRIMDATA_W, global_LENRIMDATA_W,                 &
      ! Ocean variables
     &  LENRIMO, LENRIMO_U,                                             &
      ! Variables still to be dealt with
     &  RIMFLDSA,RIMFLDSO,BOUNDFLDS,RIM_LOOKUPSA,RIM_LOOKUPSO,          &
     &  BOUND_LOOKUPSA,BOUND_LOOKUPSO,BOUND_LOOKUPSW,                   &
     &  global_LENRIMDATA_A,                                            &
     &  LENRIMDATA_A,LENRIMDATA_O
! TYPSIZE end
! TYPD1 Common block containing the ALT_N_SUBMODEL_PARTITION variables
! CALTSUBM
! TYPD1 needs access to N_SUBMODEL_PARTITION/_MAX in CSUBMODL. However,
! they are not always called in the same decks and in the right order.
! Therefore, copy the values to another file and include it from TYPD1

      INTEGER ALT_N_SUBMODEL_PARTITION

      INTEGER, PARAMETER :: ALT_N_SUBMODEL_PARTITION_MAX=4

      COMMON/CALTSUBM/ALT_N_SUBMODEL_PARTITION
! CALTSUBM end
! This file needs TYPSIZE included first

      REAL    ::  D1(LEN_TOT)       ! IN/OUT: Main data array
      LOGICAL :: LD1(LEN_TOT)       ! IN/OUT: Main data array (logical)
      INTEGER :: ID1(LEN_TOT)       ! I/OUT: Main data array (integer)

! D1_ADDR start
      ! Information for accessing D1 addressing array
      ! Number of items of info needed for each object and maximum
      ! number of objects in D1 -

      ! Number of items of information in D1 addressing array
      INTEGER,PARAMETER:: D1_LIST_LEN=17

! Names of items in D1 addressing array. Update D1_LIST_LEN above if
! items added

      ! Prognostic, Diagnostic, Secondary or other
      INTEGER,PARAMETER:: d1_object_type    = 1 ! Internal model id
      INTEGER,PARAMETER:: d1_imodl          = 2  ! Internal model id
      INTEGER,PARAMETER:: d1_section        = 3  ! Section
      INTEGER,PARAMETER:: d1_item           = 4  ! Item
      INTEGER,PARAMETER:: d1_address        = 5  ! Address in D1
      INTEGER,PARAMETER:: d1_length         = 6  ! Record length
      INTEGER,PARAMETER:: d1_grid_type      = 7  ! Grid type
      INTEGER,PARAMETER:: d1_no_levels      = 8  ! Number of levels

      ! Stash list number for diags. -1 for progs
      INTEGER,PARAMETER:: d1_stlist_no      = 9

      ! Pointer to dump header lookup table
      INTEGER,PARAMETER:: d1_lookup_ptr     = 10

      INTEGER,PARAMETER:: d1_north_code     = 11 ! Northern row
      INTEGER,PARAMETER:: d1_south_code     = 12 ! Southern row
      INTEGER,PARAMETER:: d1_east_code      = 13 ! Eastern row
      INTEGER,PARAMETER:: d1_west_code      = 14 ! Western row
      INTEGER,PARAMETER:: d1_gridpoint_code = 15 ! gridpoint info
      INTEGER,PARAMETER:: d1_proc_no_code   = 16 ! Processing Code
      INTEGER,PARAMETER:: d1_halo_type      = 17 ! Halo width type

      ! Types of items for d1_type

      INTEGER,PARAMETER:: prognostic = 0
      INTEGER,PARAMETER:: diagnostic = 1
      INTEGER,PARAMETER:: secondary  = 2
      INTEGER,PARAMETER:: other      = 3

! D1_ADDR end
      ! D1 addressing array and number of objects in each submodel
      INTEGER :: D1_ADDR(D1_LIST_LEN,N_OBJ_D1_MAX,                      &
     &  ALT_N_SUBMODEL_PARTITION)

      INTEGER :: NO_OBJ_D1(ALT_N_SUBMODEL_PARTITION_MAX)

      COMMON/common_D1_ADDRESS/ NO_OBJ_D1
! TYPD1 end
! TYPDUMA needs TYPSIZE included first
!LL  Model            Modification history
!LL version  Date
!LL   4.1    21/03/96 Added arrays to hold local lengths and addresses
!LL                   for mpp code
!LL   5.0    07/05/99 Introduce meaningful parameter names for real
!LL                   constants header. R Rawlins
!LL
!L --------------- Dump headers (atmosphere)-------------
      INTEGER :: A_FIXHD(LEN_FIXHD)    ! fixed length header
      INTEGER :: A_INTHD(A_LEN_INTHD)  ! integer header
      INTEGER :: A_CFI1(A_LEN_CFI1+1)  ! compress field index
      INTEGER :: A_CFI2(A_LEN_CFI2+1)  ! compress field index
      INTEGER :: A_CFI3(A_LEN_CFI3+1)  ! compress field index

      REAL::A_REALHD(A_LEN_REALHD)                    ! real header
      REAL::A_LEVDEPC(A_LEN1_LEVDEPC*A_LEN2_LEVDEPC+1)! level  dep const
      REAL::A_ROWDEPC(A_LEN1_ROWDEPC*A_LEN2_ROWDEPC+1)! row    dep const
      REAL::A_COLDEPC(A_LEN1_COLDEPC*A_LEN2_COLDEPC+1)! column dep const
      REAL::A_FLDDEPC(A_LEN1_FLDDEPC*A_LEN2_FLDDEPC+1)! field  dep const
      REAL::A_EXTCNST(A_LEN_EXTCNST+1)                ! extra constants
      REAL::A_DUMPHIST(LEN_DUMPHIST+1)                ! temp hist file

      ! Meaningful parameter names for integer constants header:
! ----------------------- include file: IHEADAPM -----------------------
! Description: Meaningful parameter names to index A_INTHD array in
!              atmosphere dump, ie INTEGER CONSTANTS, and reduce magic
!              numbers in code.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Original version R. Rawlins
!  5.2  07/03/01  Add height generation method P.Selwood
      INTEGER,PARAMETER:: ih_a_step          = 1  ! Timestep no.
      INTEGER,PARAMETER:: ih_rowlength       = 6  ! No. of points E-W
      INTEGER,PARAMETER:: ih_rows            = 7  ! No. of points N-S

      ! No. of model levels (0=surface)
      INTEGER,PARAMETER:: ih_model_levels    = 8

      ! No. of model levels with moisture
      INTEGER,PARAMETER:: ih_wet_levels      = 9

      ! No. of deep soil temperature levels
      INTEGER,PARAMETER:: ih_soilT_levels    = 10

      INTEGER,PARAMETER:: ih_cloud_levels    = 11 ! No. of cloud levels
      INTEGER,PARAMETER:: ih_tracer_levels   = 12 ! No. of tracer levels

      ! No. of boundary layer levels
      INTEGER,PARAMETER:: ih_boundary_levels = 13
      INTEGER,PARAMETER:: ih_N_types         = 15 ! No. of field types

       ! Height generation method
      INTEGER,PARAMETER:: ih_height_gen      = 17

      ! First rho level at which height is constant
      INTEGER,PARAMETER:: ih_1_c_rho_level   = 24

      INTEGER,PARAMETER:: ih_land_points     = 25 ! No. of land points
      INTEGER,PARAMETER:: ih_ozone_levels    = 26 ! No. of ozone levels

      ! No. of deep soil moisture levels
      INTEGER,PARAMETER:: ih_soilQ_levels    = 28

      ! Number of convective cloud levels
      INTEGER,PARAMETER:: ih_convect_levels  = 34
      INTEGER,PARAMETER:: ih_rad_step        = 35 ! Radiation timestep
      INTEGER,PARAMETER:: ih_AMIP_flag       = 36 ! Flag for AMIP run
      INTEGER,PARAMETER:: ih_AMIP_year       = 37 ! First AMIP year
      INTEGER,PARAMETER:: ih_AMIP_month      = 38 ! First AMIP month
      INTEGER,PARAMETER:: ih_AMIP_day        = 49 ! First AMIP day
      INTEGER,PARAMETER:: ih_ozone_month     = 40 ! Current ozone month
      INTEGER,PARAMETER:: ih_SH_zonal_quad   = 41 ! L_SH_zonal_quadratics
      INTEGER,PARAMETER:: ih_SH_zonal_begin  = 42 ! SH_zonal_begin
      INTEGER,PARAMETER:: ih_SH_zonal_period = 43 ! SH_zonal_period
      INTEGER,PARAMETER:: ih_SH_level_weight = 44 ! SuHe_level_weight
      INTEGER,PARAMETER:: ih_SH_sigma_cutoff = 45 ! SuHe_sigma_cutoff
      INTEGER,PARAMETER:: ih_friction_time   = 46 ! frictional_timescale

! IHEADAPM end
      ! Meaningful parameter names for real constants header:
! ----------------------- include file: RHEADAPM -----------------------
! Description: Meaningful parameter names to index A_REALHD array in
!              atmosphere dump, ie REAL CONSTANTS, and reduce magic
!              numbers in code.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Original version R. Rawlins
!  5.1  17/03/00  Change rh_tot_fluxes to rh_tot_m_init R A Stratton
!  5.1  28/04/00   Extra pressure level above top theta level
!                  p_top_theta_level removed                A.Malcolm

      ! East-West   grid spacing in degrees
      INTEGER,PARAMETER:: rh_deltaEW         = 1

      ! North-South grid spacing in degrees
      INTEGER,PARAMETER:: rh_deltaNS         = 2

      ! Latitude  of first p point in degrees
      INTEGER,PARAMETER:: rh_baselat         = 3

      ! Longitude of first p point in degrees
      INTEGER,PARAMETER:: rh_baselong        = 4

      ! Latitude  of rotated N pole in degrees
      INTEGER,PARAMETER:: rh_rotlat          = 5

      ! Longitude of rotated N pole in degrees
      INTEGER,PARAMETER:: rh_rotlong         = 6

      ! Height of top theta level (m)
      INTEGER,PARAMETER:: rh_z_top_theta     =16

      ! total moisture of the atmosphere
      INTEGER,PARAMETER:: rh_tot_m_init      =18

      ! total mass of atmosphere
      INTEGER,PARAMETER:: rh_tot_mass_init   =19

      ! total energy of atmosphere
      INTEGER,PARAMETER:: rh_tot_energy_init =20

      ! energy correction = energy drift
      INTEGER,PARAMETER:: rh_energy_corr     =21

! RHEADAPM end
      ! Meaningful parameter names for fixed header:
! ----------------------- include file: FHEADAPM -----------------------
! Description: Meaningful parameter names to index A_FIXHD array in
!              atmosphere dump, ie REAL CONSTANTS, and reduce magic
!              numbers in code.
 
      ! Start of Row Dependent Constant
      INTEGER,PARAMETER:: fh_RowDepCStart   = 115

      ! Start of Col Dependent Constant
      INTEGER,PARAMETER:: fh_ColDepCStart   = 120

! FHEADAPM end
      ! PP headers

      INTEGER :: A_LOOKUP(LEN1_LOOKUP,A_LEN2_LOOKUP) ! lookup heads
      INTEGER :: A_MPP_LOOKUP(MPP_LEN1_LOOKUP,A_LEN2_LOOKUP)
      INTEGER :: a_ixsts(len_a_ixsts)     ! stash index array

      REAL    :: a_spsts(len_a_spsts)     ! atmos stash array
! TYPDUMA end
!-----------------Start of TYPCONA------------------------------------
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field sin_u_latitude. J F Thomson.
!  5.0   18/05/99  Removed B-grid dynamics terms, replaced by C-P C-grid
!                  dynamics terms. M.L.Gallani
!  5.1   25/02/00  Add minor trig variables to allow initialisation
!                  in Setcona instead of Atm_step. R Rawlins
!  5.3   01/10/01  Add fields for chequerboard radiation. S Cusack
! 6.1  04/08/04  Add arrays  diff_coeff_u, diff_coeff_v
!                                             Terry Davies
! 6.2  04/08/04  Add true_latitude. Yongming Tang
! 6.2  14/02/06  Add arrays for variable resolution  Terry Davies
! CMAXSIZE should be included first.

      ! Constants for ATMOSPHERE model.
      !  Constants for routines independent of resolution.
! CCONSTS start
! Description:
!    This file contains declarations for derived constants within
!   the atmospheric model. Where necessary PARAMETERS are defined to
!   dimension these constants. All constants are placed in the common
!   block CDERIVED, except hardwired constants, e.g. ETA_SPLIT and LENs.
!   file CMAXSIZE must be included first
!
!   The derived constants are calculated in the routine SETCONA1.
!
! PA, WI      <- programmer of some or all of previous code or changes
!
!  Model            Modification history from model version 3.0:
! version  Date
! 3.2   26/03/93  Remove resolution dependent variables for dynamic
!                 allocation. R.Rawlins
!   4.0   20/07/95  Sizes of tables expanded for version 3A
!                   of the radiation in sections 1 and 2. J.M. Edwards
!   5.0   21/06/99  Remove obsolete constants, for C-P C-grid dynamics.
!                   M L Gallani
!   5.1   25/02/00  Replace Data ETA_SPLIT by hard-wire h_split in
!                   Readlsta. Also remove obsolete references to
!                   radiation tables. R Rawlins
! Logical component: F011

      ! No of cloud types ie low/med/high
      INTEGER, PARAMETER :: NUM_CLOUD_TYPES = 3

      ! derived constants:
      INTEGER :: LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER :: LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER :: MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER :: HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      ! height values to split model levels into l/m/h cloud
      REAL ::    h_split(NUM_CLOUD_TYPES+1)

      LOGICAL :: ELF                ! T if atmosphere model on LAM grid

      ! Constants for dynamics output independent of resolution but
      ! dependent on choice of levels for output.
      REAL :: REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)

      COMMON /CDERIVED/                                                 &
     &  h_split,LOW_BOT_LEVEL,LOW_TOP_LEVEL,MED_BOT_LEVEL,MED_TOP_LEVEL,&
     &  HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,ELF,REQ_THETA_PV_LEVS
! CCONSTS end

! typcona.h contained constants for the atmosphere.
! As of vn6.6 these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, AD_MASK_TROP_MOD, ROT_COEFF_MOD

! The following common block does not correpsond to the constants
! specified by (and subsequently commented out from) argcona.h, 
! and so it has remained here - 
! It could be moved to a more appropriate place in future

      ! Trigonometric co-ordinates in radians
      REAL:: Delta_lambda       ! EW (x) grid spacing in radians
      REAL:: Delta_phi          ! NS (y) grid spacing in radians
      REAL:: Base_phi           ! Lat of first theta point in radians
      REAL:: Base_lambda        ! Long of first theta point in radians
      REAL:: lat_rot_NP         ! Real lat of 'pseudo' N pole in radians
      REAL:: long_rot_NP        ! Real long of 'pseudo' N pole in radians

      COMMON/cderv_trig/                                                &
     &  Delta_lambda,Delta_phi,Base_phi,Base_lambda,                    &
     &  lat_rot_NP,long_rot_NP
!-------------End of TYPCONA---------------------------------------
! TYPLNDM
! Formerly integral part of TYPCONA, the variables below have been
! separated from the rest of TYPCONA as they are required by some
! of the Ocean routines in the Ocean-Atmosphere configuration of
! the UM whilest TYPCONA is not.

      ! Primary Arrays
      INTEGER::land_points     ! No. of land points  (can be 0)
      INTEGER::land_ice_points ! Number of land ice points
      INTEGER::soil_points     ! Number of soil points

!     INTEGER::land_index (land_field) ! set from land_sea_mask
!     INTEGER::land_ice_index (land_field)  ! Array of land ice points.
!     INTEGER::soil_index(land_field)       ! Array of soil points.
! sza fix conflict case when land_field=0
      INTEGER::land_index (max(1,land_field)) ! set from land_sea_mask
      INTEGER::land_ice_index (max(1,land_field))  ! Array of land ice points.
      INTEGER::soil_index(max(1,land_field))       ! Array of soil points.

      ! Gets some sizes transported around the model :
      COMMON /land_soil_dimensions/                                     &
     &  land_points , land_ice_points , soil_points

! TYPLNDA end
! CTRACERA start
!  Vn    Date    Modification History
! 6.1  23/06/04  Prognostic tracers now in section 33, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.
! 6.2  13/07/05  Also increase A_MAX_TRVARS to 150. R Barnes.
! 6.2  10/11/05  UKCA tracers put into section 34, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.

      ! First atmospheric tracer (STASH No)
      INTEGER,PARAMETER:: A_TRACER_FIRST = 1

      ! Last atmospheric tracer  (STASH No)
      INTEGER,PARAMETER:: A_TRACER_LAST = 150

      ! Maximum number of atmospheric tracers
      INTEGER,PARAMETER:: A_MAX_TRVARS  = 150
!     First UKCA tracer (STASH No)
      INTEGER,PARAMETER:: A_UKCA_FIRST = 1
!     Last UKCA tracer  (STASH No)
      INTEGER,PARAMETER:: A_UKCA_LAST = 150
!     Maximum number of UKCA tracers
      INTEGER,PARAMETER:: A_MAX_UKCAVARS  = 150

      ! Index to relative position.
      ! A_TR_INDEX(N) gives position in JTRACER for tracer number N.
      ! Set in SET_ATM_POINTERS.
      ! A_TR_INDEX(N) is the position, in the list of tracers
      ! actually present in D1, that tracer number N (in the list
      ! of all tracers selectable from the user interface) occupies,
      ! if it is present.
      ! If tracer number N is absent then A_TR_INDEX(N) is undefined.

      INTEGER :: A_TR_INDEX(A_MAX_TRVARS)
      INTEGER :: A_UKCA_INDEX(A_MAX_UKCAVARS)   ! Similar to A_TR_INDEX

      COMMON/ATRACER/A_TR_INDEX,A_UKCA_INDEX

! CTRACERA end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Start typ_atm_fields.h

! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "typptra.h", and requires that "ctracera.h" is
!  also included.  By 7.1 these should all be their natural shape instead
!  of all being 1D arrays.
!
! Current Code Owner: A. Treshansky
!

      ! 1.1: Data variables stored in primary space.
      real :: U(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: V(1-offx:row_length+offx,1-offy:n_rows+offy,model_levels)
      real :: W(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)
      real :: RHO(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: THETA(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      real :: Q(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: QCL(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: QCF(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: QCF2(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: QRAIN(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: QGRAUP(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)

! EXNER_RHO_LEVELS is still 1D; change for 7.1
      ! Exner pressure on rho levels
      real :: EXNER_RHO_LEVELS(((2*offx)+row_length)*((2*offy)+rows)*(model_levels+1))   

      real :: U_ADV(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: V_ADV(1-offx:row_length+offx,1-offy:n_rows+offy,model_levels)
      real :: W_ADV(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

! P is still 1D; change for 7.1
      ! 1.2: Data variables stored in secondary space.
      real :: P(((2*offx)+row_length)*((2*offy)+rows)*(model_levels+1))

      ! Pressure on theta levels
      real :: P_THETA_LEVELS(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      ! Exner pressure on theta levels
      real :: EXNER_THETA_LEVELS(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      ! 1.3: Cloud Fields
      real :: CCW_RAD(row_length,rows,wet_levels)
      ! CCA's size is dependant on L_3D_CCA
      ! N_CCA_LEV will be set to the correct value (either wet_levels or 1)
      real :: CCA(row_length,rows,n_cca_lev)
      real :: CF_AREA(row_length,rows,wet_levels)
      real :: CF_BULK(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: CF_LIQUID(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: CF_FROZEN(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)

      ! 1.4: Soil Ancillary fields
      real :: DEEP_SOIL_TEMP(land_field,sm_levels)
      real :: SMCL(land_field,sm_levels)
      real :: STHU(land_field,sm_levels)
      real :: STHF(land_field,sm_levels)

      ! 1.5: Radiation Increments
      real :: SW_INCS(row_length,rows,0:model_levels+1)  ! SW radiation increments
      real :: LW_INCS(row_length,rows,0:model_levels)    ! LW radiation increments
! PAR radiation increment
      real :: DIRPAR(row_length,rows)

! OZONE fields are still 1D; change for 7.1
      ! 1.6: Ozone
      real :: O3(row_length*rows*ozone_levels)
!  tropopause-based ozone
      real :: TPPSOZONE(row_length*rows*ozone_levels)
!  Cariolle ozone tracer variables
      real :: OZONE_TRACER(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: O3_PROD_LOSS(1,rows,model_levels)
      real :: O3_P_L_VMR(1,rows,model_levels)
      real :: O3_VMR(1,rows,model_levels)
      real :: O3_P_L_TEMP(1,rows,model_levels)
      real :: O3_TEMP(1,rows,model_levels)
      real :: O3_P_L_COLO3(1,rows,model_levels)
      real :: O3_COLO3(1,rows,model_levels)
      
! TRACERS are still 1D; change for 7.1
      ! 1.7: Tracer and aerosol fields
      ! TRACERS are dealt w/ differently
      ! these are the maximum sizes:
      real :: TRACER(tr_levels*theta_off_size*(a_tracer_last-a_tracer_first))
      real :: TRACER_UKCA(tr_levels*theta_off_size*(a_ukca_last-a_ukca_first))
      real :: TRACER_LBC(LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_vars)
      real :: TRACER_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_vars)

      real :: MURK_SOURCE(row_length,rows,model_levels)
      real :: MURK(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV4(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV5(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV6(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SO2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DMS(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SO4_AITKEN(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SO4_ACCU(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SO4_DISS(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: H2O2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: NH3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SOOT_NEW(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SOOT_AGD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SOOT_CLD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: BMASS_NEW(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: BMASS_AGD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: BMASS_CLD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: OCFF_NEW(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: OCFF_AGD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: OCFF_CLD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SO2_NATEM(row_length,rows,model_levels)
      real :: OH(row_length,rows,model_levels)
      real :: HO2(row_length,rows,model_levels)
      real :: H2O2_LIMIT(row_length,rows,model_levels)
      real :: O3_CHEM(row_length,rows,model_levels)
      real :: CO2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: CH4_STOCH(row_length,rows,model_levels)
      real :: O3_STOCH(row_length,rows,model_levels)

! USER_MULT<N> are still 1D; change for 7.1
! 1.8: Multi-level user ancillary fields
      real :: USER_MULT1(row_length*rows*model_levels)
      real :: USER_MULT2(row_length*rows*model_levels)
      real :: USER_MULT3(row_length*rows*model_levels)
      real :: USER_MULT4(row_length*rows*model_levels)
      real :: USER_MULT5(row_length*rows*model_levels)
      real :: USER_MULT6(row_length*rows*model_levels)
      real :: USER_MULT7(row_length*rows*model_levels)
      real :: USER_MULT8(row_length*rows*model_levels)
      real :: USER_MULT9(row_length*rows*model_levels)
      real :: USER_MULT10(row_length*rows*model_levels)
      real :: USER_MULT11(row_length*rows*model_levels)
      real :: USER_MULT12(row_length*rows*model_levels)
      real :: USER_MULT13(row_length*rows*model_levels)
      real :: USER_MULT14(row_length*rows*model_levels)
      real :: USER_MULT15(row_length*rows*model_levels)
      real :: USER_MULT16(row_length*rows*model_levels)
      real :: USER_MULT17(row_length*rows*model_levels)
      real :: USER_MULT18(row_length*rows*model_levels)
      real :: USER_MULT19(row_length*rows*model_levels)
      real :: USER_MULT20(row_length*rows*model_levels)

!     CABLE
      real :: TSOIL_TILE(land_field,ntiles,st_levels)
      real :: SMCL_TILE(land_field,ntiles,sm_levels)
! Not used MRD
!      real :: STHU_TILE(land_field,ntiles,sm_levels)
      real :: STHF_TILE(land_field,ntiles,sm_levels)
      real :: SNOW_DEPTH3L(land_field,ntiles,3)
      real :: SNOW_MASS3L(land_field,ntiles,3)
      real :: SNOW_TMP3L(land_field,ntiles,3)
      real :: SNOW_RHO3L(land_field,ntiles,3)
      real :: SNOW_RHO1L(land_field,ntiles)
      real :: SNOW_AGE(land_field,ntiles)
      real :: SNOW_FLG3L(land_field,ntiles)

! Lestevens Sept 2012 - For CASA-CNP
      real :: CPOOL_TILE(land_field,ntiles,10)
      real :: NPOOL_TILE(land_field,ntiles,10)
      real :: PPOOL_TILE(land_field,ntiles,12)
      real :: SOIL_ORDER(land_field)
      real :: NIDEP(land_field)
      real :: NIFIX(land_field)
      real :: PWEA(land_field)
      real :: PDUST(land_field)
      real :: GLAI(land_field,ntiles)
      real :: PHENPHASE(land_field,ntiles)

      real :: cable_lai(land_field,ntiles)
      
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
      real :: OROG_LBC(LENRIMA(fld_type_p,halo_type_extended,1))
      real :: U_LBC(LENRIMA(fld_type_u,halo_type_extended,1),model_levels)
      real :: V_LBC(LENRIMA(fld_type_v,halo_type_extended,1),model_levels)
      real :: W_LBC(LENRIMA(fld_type_p,halo_type_extended,1),0:model_levels)
      real :: RHO_LBC(LENRIMA(fld_type_p,halo_type_extended,1),model_levels)
      real :: THETA_LBC(LENRIMA(fld_type_p,halo_type_extended,1),model_levels)
      real :: Q_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCL_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCF_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCF2_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QRAIN_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QGRAUP_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_BULK_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_LIQUID_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_FROZEN_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: EXNER_LBC(LENRIMA(fld_type_p,halo_type_extended,1),model_levels+1)
      real :: U_ADV_LBC(LENRIMA(fld_type_u,halo_type_extended,1),model_levels)
      real :: V_ADV_LBC(LENRIMA(fld_type_v,halo_type_extended,1),model_levels)
      real :: W_ADV_LBC(LENRIMA(fld_type_p,halo_type_extended,1),0:model_levels)
      real :: MURK_LBC(LENRIMA(fld_type_p,halo_type_single,1),model_levels)
      
      ! Lateral Boundary Condition tendencies
      real :: U_LBC_TEND(LENRIMA(fld_type_u,halo_type_extended,1),model_levels)
      real :: V_LBC_TEND(LENRIMA(fld_type_v,halo_type_extended,1),model_levels)
      real :: W_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),0:model_levels)
      real :: RHO_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),model_levels)
      real :: THETA_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),model_levels)
      real :: Q_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCL_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCF_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCF2_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QRAIN_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QGRAUP_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_BULK_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_LIQUID_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_FROZEN_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: EXNER_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),model_levels+1)
      real :: U_ADV_LBC_TEND(LENRIMA(fld_type_u,halo_type_extended,1),model_levels)
      real :: V_ADV_LBC_TEND(LENRIMA(fld_type_v,halo_type_extended,1),model_levels)
      real :: W_ADV_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),0:model_levels)
      real :: MURK_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),model_levels)

      ! 2: Scalar Variables

      ! 2.1: Data variables stored in primary space.
      real :: TSTAR(row_length,rows)
      real :: LAND(row_length,rows)
      real :: TSTAR_ANOM(row_length,rows)
      ! 2.15: Fields for coastal tiling
      real :: FRAC_LAND(land_field)
      real :: TSTAR_LAND(row_length,rows)
      real :: TSTAR_SEA(row_length,rows)
      real :: TSTAR_SICE(row_length,rows)
      ! Set pointers for sea-ice and land albedos
      real :: SICE_ALB(row_length,rows)
      real :: LAND_ALB(row_length,rows)

      ! 2.2: Data variables stored in secondary space.

      real :: PSTAR(row_length,rows)

      ! 2.3: Cloud fields
      real :: LCBASE(row_length,rows)
      real :: CCB(row_length,rows)
      real :: CCT(row_length,rows)
      real :: CCLWP(row_length,rows)

      ! 2.4: Boundary layer fields

      real :: ZH(row_length,rows)

      ! Standard deviation of turbulent fluctuations of layer 1 temperature
      real :: T1_SD(row_length,rows)

      ! Standard deviation of turbulent fluctuations of layer 1 humidity
      real :: Q1_SD(row_length,rows)

      ! Number of model levels in the  turbulently mixed layer
      real :: NTML(row_length,rows)

      ! Top level for turb mixing in any decoupled Sc layer
      real :: NTDSC(row_length,rows)

      ! Bottom level for turb mixing in any decoupled Sc layer
      real :: NBDSC(row_length,rows)

      real :: CUMULUS(row_length,rows)

      ! 2.4: Soil Ancillary fields

      real :: SAT_SOILW_SUCTION(land_field)
      real :: THERM_CAP(land_field)
      real :: THERM_COND(land_field)
      real :: VOL_SMC_CRIT(land_field)
      real :: VOL_SMC_WILT(land_field)
      real :: VOL_SMC_SAT(land_field)
      real :: SAT_SOIL_COND(land_field)
      real :: CLAPP_HORN(land_field)

      ! 2.5: Vegetation Ancillary fields

      real :: Z0(row_length,rows)
      real :: CANOPY_WATER(land_field)
      real :: SURF_CAP(row_length,rows)
      real :: SURF_RESIST(land_field)
      real :: ROOT_DEPTH(land_field)
      real :: INFILT(land_field)
      real :: VEG_FRAC(land_field)
      real :: LAI(land_field)
      real :: CANHT(land_field)
      real :: SFA(land_field)
      real :: MDSA(land_field)
      real :: GS(land_field)

      ! 2.6: Orographic Ancillary fields

      real :: OROGRAPHY(row_length,rows)
      real :: OROG_SD(land_field)
      real :: OROG_SIL(land_field)
      real :: OROG_HO2(land_field)
      real :: OROG_GRAD_X(land_field)
      real :: OROG_GRAD_Y(land_field)
      real :: OROG_GRAD_XX(land_field)
      real :: OROG_GRAD_XY(land_field)
      real :: OROG_GRAD_YY(land_field)

      ! 2.7: Sea/Sea Ice

      real :: U_SEA(row_length,rows)
      real :: V_SEA(row_length,n_rows)
      real :: U_0_P(row_length,rows)
      real :: V_0_P(row_length,rows)

! these fields are targets b/c there are pointers w/in ATM_STEP
! which point to one or another of them depending on the ice category selection
      real, target :: TI(row_length,rows,1)
      real, target :: ICE_FRACTION(row_length,rows,1)
      real, target :: ICE_THICKNESS(row_length,rows,1)
      real, target :: TI_CAT(row_length,rows,nice) 
      real, target :: ICE_FRACT_CAT(row_length,rows,nice)
      real, target :: ICE_THICK_CAT(row_length,rows,nice)

      ! 2.8: Snow

      real :: SNODEP(row_length,rows)
      real :: SNODEP_SEA(row_length,rows)
      real :: SNODEP_SEA_CAT(row_length,rows,nice)
      real :: CATCH_SNOW(land_field)
      real :: SNOW_GRND(land_field)
      ! SNSOOT may not be used as of vn6.6
      real :: SNSOOT(row_length,rows)  ! Snow soot content

      ! 2.9: aerosol emission fields, including mineral dust parent soil props

      real :: SOIL_CLAY(row_length,rows)
      real :: SOIL_SILT(row_length,rows)
      real :: SOIL_SAND(row_length,rows)
      real :: DUST_MREL1(row_length,rows)
      real :: DUST_MREL2(row_length,rows)
      real :: DUST_MREL3(row_length,rows)
      real :: DUST_MREL4(row_length,rows)
      real :: DUST_MREL5(row_length,rows)
      real :: DUST_MREL6(row_length,rows)

      real :: SO2_EM(row_length,rows)
      real :: DMS_EM(row_length,rows)
      real :: SO2_HILEM(row_length,rows)
      real :: NH3_EM(row_length,rows)
      real :: SOOT_EM(row_length,rows)
      real :: SOOT_HILEM(row_length,rows)
      real :: BMASS_EM(row_length,rows)
      real :: BMASS_HILEM(row_length,rows)
      real :: OCFF_EM(row_length,rows)
      real :: OCFF_HILEM(row_length,rows)
      real :: DMS_CONC(row_length,rows)
      real :: DMS_OFLUX(row_length,rows)

! USER_ANC<N> fields are still 1D; change for 7.1
      ! 2.10: User ancillary fields
      real :: USER_ANC1(row_length*rows*1)
      real :: USER_ANC2(row_length*rows*1)
      real :: USER_ANC3(row_length*rows*1)
      real :: USER_ANC4(row_length*rows*1)
      real :: USER_ANC5(row_length*rows*1)
      real :: USER_ANC6(row_length*rows*1)
      real :: USER_ANC7(row_length*rows*1)
      real :: USER_ANC8(row_length*rows*1)
      real :: USER_ANC9(row_length*rows*1)
      real :: USER_ANC10(row_length*rows*1)
      real :: USER_ANC11(row_length*rows*1)
      real :: USER_ANC12(row_length*rows*1)
      real :: USER_ANC13(row_length*rows*1)
      real :: USER_ANC14(row_length*rows*1)
      real :: USER_ANC15(row_length*rows*1)
      real :: USER_ANC16(row_length*rows*1)
      real :: USER_ANC17(row_length*rows*1)
      real :: USER_ANC18(row_length*rows*1)
      real :: USER_ANC19(row_length*rows*1)
      real :: USER_ANC20(row_length*rows*1)

      ! Tracer fluxes - kdcorbin, 05/10
      real :: TRACER_FLUX1(row_length,rows)
      real :: TRACER_FLUX2(row_length,rows)
      real :: TRACER_FLUX3(row_length,rows)
      real :: TRACER_FLUX4(row_length,rows)
      real :: TRACER_FLUX5(row_length,rows)
      real :: TRACER_FLUX6(row_length,rows)
      real :: TRACER_FLUX7(row_length,rows)
      real :: TRACER_FLUX8(row_length,rows)
      real :: TRACER_FLUX9(row_length,rows)
      real :: TRACER_FLUX10(row_length,rows)
      real :: TRACER_FLUX11(row_length,rows)
      real :: TRACER_FLUX12(row_length,rows)
      real :: TRACER_FLUX13(row_length,rows)
      real :: TRACER_FLUX14(row_length,rows)
      real :: TRACER_FLUX15(row_length,rows)
      real :: TRACER_FLUX16(row_length,rows)
      real :: TRACER_FLUX17(row_length,rows)
      real :: TRACER_FLUX18(row_length,rows)
      real :: TRACER_FLUX19(row_length,rows)
      real :: TRACER_FLUX20(row_length,rows)

      !   2.11: Store arrays for energy correction calculation
      real :: NET_FLUX(row_length,rows)
      real :: NET_MFLUX(row_length,rows)

      !   2.12: Tiled Vegetation and Triffid fields
      real :: FRAC_TYP(land_field)
      real :: FRAC_CON1(land_field)  ! fraction of broadleaf tree
      real :: FRAC_CON2(land_field)  ! fraction of needleleaf tree
      real :: FRAC_CON3(land_field)  ! fraction of C3 grass
      real :: FRAC_CON4(land_field)  ! fraction of C4 grass
      real :: FRAC_CON5(land_field)  ! fraction of shrub
      real :: FRAC_CON6(land_field)  ! fraction of urban
      real :: FRAC_CON7(land_field)  ! fraction of water
      real :: FRAC_CON8(land_field)  ! fraction of soil
      real :: FRAC_CON9(land_field)  ! fraction of ice
      real :: LAI_PFT(land_field)
      real :: CANHT_PFT(land_field)
      real :: DISTURB_VEG(land_field)
      real :: SOIL_ALB(land_field)
      real, target :: SOIL_CARB(land_field)
      real, target :: SOIL_CARB1(land_field)
      real, target :: SOIL_CARB2(land_field)
      real, target :: SOIL_CARB3(land_field)
      real, target :: SOIL_CARB4(land_field)
      real :: NPP_PFT_ACC(land_field)
      real :: G_LF_PFT_ACC(land_field)
      real :: G_PHLF_PFT_ACC(land_field)
      real :: RSP_W_PFT_ACC(land_field)
      real, target :: RSP_S_ACC(land_field)
      real, target :: RSP_S_ACC1(land_field)
      real, target :: RSP_S_ACC2(land_field)
      real, target :: RSP_S_ACC3(land_field)
      real, target :: RSP_S_ACC4(land_field)
      real :: CAN_WATER_TILE(land_field)
      real :: CATCH_TILE(land_field)
      real :: INFIL_TILE(land_field)
      real :: RGRAIN_TILE(land_field)
      real :: SNODEP_TILE(land_field)
      real :: TSTAR_TILE(land_field) 
      real :: Z0_TILE(land_field)  
      real :: DOLR_FIELD(row_length,rows) 
      real :: LW_DOWN(row_length,rows) 
      real :: SW_TILE_RTS(land_field)

!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!      real :: TSLAB(row_length*rows*1)
!      real :: TCLIM(row_length*rows*1)
!      real :: HCLIM(row_length*rows*1)
!      real :: CHEAT(row_length*rows*1)
!      real :: OIFLX(row_length*rows*1)
!      real :: UICE(row_length*rows*1)
!      real :: VICE(row_length*n_rows*1)
!      real :: SIG11NE(row_length*rows*1)
!      real :: SIG11SE(row_length*rows*1) 
!      real :: SIG11SW(row_length*rows*1)
!      real :: SIG11NW(row_length*rows*1)
!      real :: SIG12NE(row_length*rows*1)
!      real :: SIG12SE(row_length*rows*1)
!      real :: SIG12NW(row_length*rows*1)
!      real :: SIG22NE(row_length*rows*1)
!      real :: SIG22SE(row_length*rows*1)
!      real :: SIG22SW(row_length*rows*1)
!      real :: SIG22NW(row_length*rows*1)

!   2.14: Carbon cycle fields
      real :: CO2FLUX(row_length,rows)
      real :: CO2_EMITS(row_length,rows)

!   2.15: Fields carried forward from previous version.
!         May not be required
!      real, pointer :: SURF_RESIST_NIT(:)  ! Surface resistance on
!                                    ! non-ice tiles
!      real, pointer :: ROOT_DEPTH_NIT(:)   ! Root depth on non-ice tiles
!      real, pointer :: Z0V_TYP(:)          ! Vegetative Roughness length on
!                                    ! tiles
!      real, pointer :: ICE_EDGE(:)
!      real, pointer :: OROG_TENDENCY(:)    ! Orographic tendencies
!      real, pointer :: OROG_SD_TENDENCY(:) ! Orographic variable tendency

      ! Pointers for ATMOSPHERE model constants. Scalars only.
      ! Addresses in level dependent constants array.
!      real, pointer :: ETATHETA(:)
!      real, pointer :: ETARHO(:)
!      real, pointer :: RHCRIT(:)
!      real, pointer :: SOIL_THICKNESS(:)
! these fields are still 1D; change for 7.1
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      real :: zseak_theta(0:model_levels)
      real :: Ck_theta(0:model_levels)
      real :: zseak_rho(model_levels)
      real :: Ck_rho(model_levels)   

      ! 2.16: Fields for large-scale hydrology scheme.
      real :: TI_MEAN(land_field)
      real :: TI_SIG(land_field)
      real :: FEXP(land_field)
      real :: GAMMA_INT(land_field)
      real :: WATER_TABLE(land_field)
      real :: FSFC_SAT(land_field)
      real :: F_WETLAND(land_field)

      real :: STHZW(land_field)
      real :: A_FSAT(land_field)
      real :: C_FSAT(land_field)
      real :: A_FWET(land_field)
      real :: C_FWET(land_field)

      ! 2.17: Fields for River routing.
      real :: RIV_SEQUENCE(river_row_length,river_rows)
      real :: RIV_DIRECTION(river_row_length,river_rows)
      real :: RIV_STORAGE(river_row_length,river_rows)
      real :: TOT_SURFROFF(river_row_length,river_rows)
      real :: TOT_SUBROFF(river_row_length,river_rows)
      real :: RIV_INLANDATM(land_field)
      ! Fields for grid-to-grid river routing (river routing 2A)
      real :: RIV_IAREA(row_length,rows)      ! Drainage area
      real :: RIV_SLOPE(row_length,rows)      ! Grid-cell slope
      real :: RIV_FLOWOBS1(row_length,rows)   ! Initial values of flow
      real :: RIV_INEXT(row_length,rows)      ! Flow direction (x)
      real :: RIV_JNEXT(row_length,rows)      ! Flow direction (y)
      real :: RIV_LAND(row_length,rows)       ! Land-type (land/river/sea)
      real :: RIV_SUBSTORE(row_length,rows)   ! Subsurface storage
      real :: RIV_SURFSTORE(row_length,rows)  ! Surface storage
      real :: RIV_FLOWIN(row_length,rows)     ! Surface inflow
      real :: RIV_BFLOWIN(row_length,rows)    ! Subsurface inflow

! Fields used when coupling using OASIS.
      real :: C_SOLAR(row_length,rows)         ! CPL solar radn
      real :: C_BLUE(row_length,rows)          ! CPL blue radn
      real :: C_DOWN(row_length,rows)          ! CPL downward radn
      real :: C_LONGWAVE(row_length,rows)      ! CPL lw radn
      real :: C_TAUX(row_length,rows)          ! CPL taux 
      real :: C_TAUY(row_length,rows)          ! CPL tauy 
      real :: C_WINDMIX(row_length,rows)       ! CPL WME   
      real :: C_SENSIBLE(row_length,rows)      ! CPL sensible ht flx
      real :: C_SUBLIM(row_length,rows)        ! CPL sublim rate
      real :: C_EVAP(row_length,rows)          ! CPL Evap rate
      real :: C_BOTMELTN(row_length,rows,nice) ! CPL Multi-cat bmlt
      real :: C_TOPMELTN(row_length,rows,nice) ! CPL Multi-cat tmlt
      real :: C_LSRAIN(row_length,rows)        ! CPL Lg scl rain rate
      real :: C_LSSNOW(row_length,rows)        ! CPL Lg scl snow rate
      real :: C_CVRAIN(row_length,rows)        ! CPL Cnvctv rain rate
      real :: C_CVSNOW(row_length,rows)        ! CPL Cnvctv snur rate
      real :: C_RIVEROUT(row_length,rows)      ! CPL Riv outflow->ocn                     

      real :: C_PRESS(row_length,rows)         ! CPL Surf pressure->ocn                     
      real :: C_U10(row_length,rows)         ! CPL 10m wind speed
      real :: C_V10(row_length,rows)         ! CPL 10m wind speed


! UKCA fields are still 1D; change for 7.1
      ! UKCA oxidant fields
      real :: OH_UKCA(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: HO2_UKCA(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: H2O2_UKCA(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: O3_UKCA(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      ! Aerosol climatologies
      real :: ARCLBIOG_BG(row_length,rows,model_levels)
      real :: ARCLBIOM_FR(row_length,rows,model_levels)
      real :: ARCLBIOM_AG(row_length,rows,model_levels)
      real :: ARCLBIOM_IC(row_length,rows,model_levels)
      real :: ARCLBLCK_FR(row_length,rows,model_levels)
      real :: ARCLBLCK_AG(row_length,rows,model_levels)
      real :: ARCLSSLT_FI(row_length,rows,model_levels)
      real :: ARCLSSLT_JT(row_length,rows,model_levels)
      real :: ARCLSULP_AC(row_length,rows,model_levels)
      real :: ARCLSULP_AK(row_length,rows,model_levels)
      real :: ARCLSULP_DI(row_length,rows,model_levels)
      real :: ARCLDUST_B1(row_length,rows,model_levels)
      real :: ARCLDUST_B2(row_length,rows,model_levels)
      real :: ARCLDUST_B3(row_length,rows,model_levels)
      real :: ARCLDUST_B4(row_length,rows,model_levels)
      real :: ARCLDUST_B5(row_length,rows,model_levels)
      real :: ARCLDUST_B6(row_length,rows,model_levels)
      real :: ARCLOCFF_FR(row_length,rows,model_levels)
      real :: ARCLOCFF_AG(row_length,rows,model_levels)
      real :: ARCLOCFF_IC(row_length,rows,model_levels)
      real :: ARCLDLTA_DL(row_length,rows,model_levels)

! End typ_atm_fields.h
! TYPBND - needs TYPSIZE included first
!LL
!LL  4.5  04/08/97 Add O_BDY_STEP_PREV for ocean boundary routines
!LL                Delete FLOOR_STEPSO.          C.G. Jones
!LL  5.0  28/04/99 Remove references to FLOOR variables        P.Burton
!LL  5.2  20/09/00 Removed old LBC stuff
!LL                Removed 2nd dimension (for lower boundary)
!LL                Added LOOKUP_COMP_BOUNDA             P.Burton
!LL  5.3  15/10/01 Include additional ocean boundary data    M J Bell
!LL  5.5  17/02/03 Include Wave model boundary data. D.Holmes-Bell
!LL  6.2  23/11/05  Removed all references to the wavemodel.
!LL                 T.Edwards
! CBOUND start

!  History:
!  Date      Vn     Modification
!  31/10/01  5.3    Remove RIMWEIGHTS_OROG. D. Robinson
!  17/09/02  5.4    Variables for controlling 2nd bndy file. A Clayton
!  17/02/03  5.5    Allow Wave model to use boundary code. D.Holmes-Bell
!  20/01/06  6.2    Add Current_LBC_Step. Dave Robinson
!  01/03/06  6.2    Remove RIMWEIGHTSW. Dave Robinson

      ! These 3 arrays are set by namelist read in IN_BOUND hence
      ! cannot be in argument list nor in COMMON if array lengths are
      ! passed variables. Only way seems to be to set MAX allowed
      ! lengths consistent with User Interface so that can be in COMMON
      INTEGER, PARAMETER :: MAX_BND_FLDS=4
      INTEGER, PARAMETER :: MAX_RIMWIDTH=10
      INTEGER :: BOUND_FIELDCODE(MAX_BND_FLDS)  ! Set by NAMELIST
      REAL :: RIMWEIGHTSA(MAX_RIMWIDTH)      ! Set by NAMELIST
      REAL :: RIMWEIGHTSO(MAX_RIMWIDTH)      ! Set by NAMELIST

      ! Variable for controlling LBC updating
      ! - initialised in INBOUNDA and updated in BOUNDVAL

      Integer  :: Current_LBC_Step     ! Timestep at which LBCs were
                                       ! last updated.

      ! Variables for controlling 2nd atmos boundary file. All are
      ! calculated within the code from other data.
      INTEGER :: ALBC_num              ! Number of atmos boundary file
                                       ! currently in use
      INTEGER :: ALBC2_StartTime_steps ! VT of first block of data in
                                       ! 2nd atmos boundary file, in
                                       ! steps from start of run
      INTEGER :: ALBC_SwapStep         ! Step on which to swap to 2nd
                                       ! atmos boundary file

      COMMON /BOUND_CT/ BOUND_FIELDCODE,                                &
     &                  RIMWEIGHTSA, RIMWEIGHTSO,                       &
     &                  ALBC_num, ALBC2_StartTime_steps, ALBC_SwapStep, &
     &                  Current_LBC_Step

! CBOUND end

      !  Control data calculated from namelist
      INTEGER :: RIM_STEPSA      ! Set by IN_BOUND from BOUND_FIELDCODE
      INTEGER :: RIM_STEPSO      ! Set by IN_BOUND from BOUND_FIELDCODE
      INTEGER::RIM_STEPSW      ! Set by IN_BOUND from BOUND_FIELDCODE
      INTEGER :: NBOUND_LOOKUP(2)
      INTEGER :: O_BDY_STEP_NEXT ! timestep for which next boundary data
                          ! is valid. Calculated in INBOUND / UPBOUND
      INTEGER::W_BDY_STEP_NEXT
      COMMON/CBND/                                                      &
     & RIM_STEPSA,RIM_STEPSO,O_BDY_STEP_NEXT                            &
     &  ,W_BDY_STEP_NEXT,RIM_STEPSW
! TYPBND end
! TYPSTS starts
! CSUBMODL must be included before this file
!Applicable to all configurations (except MOS variables)
!STASH related variables for describing output requests and space
!management.
!LL
!LL   AUTHOR            Rick Rawlins
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL   3.2             Code creation for Dynamic allocation
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL   3.5  Apr. 95   Sub-Models project.
!LL                  Dimensioning of various STASH arrays altered in
!LL                  accordance with internal model separation scheme.
!LL                  Arrays PPXREF, INDEX_PPXREF deleted as they are no
!LL                  longer required.
!LL                  S.J.Swarbrick
!LL
!
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS

      INTEGER :: MOS_MASK_LEN         ! Size of bit mask for MOS

      COMMON/DSIZE_AO/  MOS_MASK_LEN

! TYPSTSZ end
!LL  Comdeck: CPPXREF --------------------------------------------------
!LL
!LL  Purpose: Holds PARAMETER definitions to describe the structure of
!LL           each STASHmaster file record plus some valid entries.
!LL
!LL  Author    Dr T Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                  Add a PPXREF record for model number.
!LL  4.0   26/07/95  T.Johns.  Add codes for real/int/log data types.
!LL  3.5   10/3/94   Sub-Models project:
!LL                 List of PPXREF addressing codes augmented, in order
!LL                 to include all of the pre_STASH master information
!LL                 in the new PPXREF file.
!LL                 PPXREF_CODELEN increased to 38.
!LL                 PPXREF_IDLEN deleted - no longer relevant.
!LL                   S.J.Swarbrick
!LL  4.1   June 96  Wave model parameters included.
!LL                 ppx_ address parameters adjusted to allow for
!LL                  reading option code as 4x5 digit groups.
!LL                   S.J.Swarbrick
!LL  5.0   29/06/99  Add halo type parameter for new dynamics.
!LL                  New grid codes for LAM boundary conditions
!LL                  D.M. Goddard
!LL  5.1   07/03/00  Fixed/Free format conversion
!LL  5.2   19/09/00  Added ppx_atm_lbc_orog descriptor   P.Burton
!LL  5.3   21/08/01  Added ocean lbc descriptors.   M. J. Bell
!LL  5.3   23/07/01  Add valid pp_lbvc codes referenced in UM. R Rawlins
!LL  5.5   30/01/03  Option code increase from 20 to 30 digits thus
!LL                  requiring option code address range increase by
!LL                  2 so all subsequent addressing codes need to be
!LL                  increased by 2 to make a gap.
!LL                  W Roseblade
!LL
!LL  Logical components covered: C40
!LL
!-----------------------------------------------------------------------
! Primary file record definition
      ! length of ID in a record
      Integer, Parameter :: PPXREF_IDLEN      = 2

      ! total length of characters *WARNING* must be multiple of 4
      ! to avoid overwriting
      Integer, Parameter :: PPXREF_CHARLEN    = 36

      ! number of packing profiles
      Integer, Parameter :: PPXREF_PACK_PROFS = 10

      ! total length of codes = no. of codes (excluding profs)
      ! + pack_profs
      Integer, Parameter :: PPXREF_CODELEN    = 33 + PPXREF_PACK_PROFS

! Derived file record sizes
      ! Assume that an integer is at least 4 bytes long. Wastes some
      ! space on an 8 byte machine.
      ! ppx_charword = 9.
      Integer, Parameter :: PPX_CHARWORD      = ((PPXREF_CHARLEN+3)/4)

      ! read buffer record length
      Integer, Parameter :: PPX_RECORDLEN = PPX_CHARWORD+PPXREF_CODELEN
!
!-----------------------------------------------------------------------
! Addressing codes within PPXREF
      Integer, Parameter ::  ppx_model_number   = 1  ! Model number
                                                     ! address
      Integer, Parameter ::  ppx_section_number = 2  ! Section number
                                                     ! address
      Integer, Parameter ::  ppx_item_number    = 3  ! Item number
                                                     ! address
      Integer, Parameter ::  ppx_version_mask   = 4  ! Version mask
                                                     ! address
      Integer, Parameter ::  ppx_space_code     = 5  ! Space code
                                                     ! address
      Integer, Parameter ::  ppx_timavail_code  = 6  ! Time availability
                                                     !  code  address
      Integer, Parameter ::  ppx_grid_type      = 7  ! Grid type code
                                                     ! address
      Integer, Parameter ::  ppx_lv_code        = 8  ! Level type code
                                                     ! address
      Integer, Parameter ::  ppx_lb_code        = 9  ! First level code
                                                     !  address
      Integer, Parameter ::  ppx_lt_code        =10  ! Last level code
                                                     ! address
      Integer, Parameter ::  ppx_lev_flag       =11  ! Level compression
                                                     !  flag  address
      Integer, Parameter ::  ppx_opt_code       =12  ! Sectional option
                                                     ! code  address
      Integer, Parameter ::  ppx_pt_code        =18  ! Pseudo dimension
                                                     ! type  address
      Integer, Parameter ::  ppx_pf_code        =19  ! First pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_pl_code        =20  ! Last pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_ptr_code       =21  ! Section 0 point-
                                                     ! back code address
      Integer, Parameter ::  ppx_dump_packing   =22  ! Dump packing code
                                                     ! address
      Integer, Parameter ::  ppx_lbvc_code      =23  ! PP LBVC code
                                                     ! address
      Integer, Parameter ::  ppx_rotate_code    =24  ! Rotation code
                                                     ! address
      Integer, Parameter ::  ppx_field_code     =25  ! PP field code
                                                     ! address
      Integer, Parameter ::  ppx_user_code      =26  ! User code address
      Integer, Parameter ::  ppx_meto8_levelcode=27  ! CF level code
                                                     ! address
      Integer, Parameter ::  ppx_meto8_fieldcode=28  ! CF field code
                                                     ! address
      Integer, Parameter ::  ppx_cf_levelcode   =27
      Integer, Parameter ::  ppx_cf_fieldcode   =28
      Integer, Parameter ::  ppx_base_level     =29  ! Base level code
                                                     ! address
      Integer, Parameter ::  ppx_top_level      =30  ! Top level code
                                                     ! address
      Integer, Parameter ::  ppx_ref_lbvc_code  =31  ! Ref level LBVC
                                                     ! code address
      Integer, Parameter ::  ppx_data_type      =32  ! Data type code
                                                     ! address
      Integer, Parameter ::  ppx_halo_type      =33
      Integer, Parameter ::  ppx_packing_acc    =34  ! Packing accuracy
                                                     ! code  address
      Integer, Parameter ::  ppx_pack_acc       =34  ! Must be last:


                                                 ! multiple pack_acc to
                                                 ! fill up remaining
                                                 ! array elements


!-------------------------------------------------------------------
! Valid grid type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_atm_nonstd=0      ! Non-standard atmos
                                                  ! grid
      Integer, Parameter :: ppx_atm_tall=1        ! All T points (atmos)
      Integer, Parameter :: ppx_atm_tland=2       ! Land-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tsea=3        ! Sea-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tzonal=4      ! Zonal field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_tmerid=5      ! Merid field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_uall=11       ! All u points (atmos)
      Integer, Parameter :: ppx_atm_uland=12      ! Land-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_usea=13       ! Sea-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_uzonal=14     ! Zonal field at u
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_umerid=15     ! Merid field at u
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_scalar=17     ! Scalar (atmos)
      Integer, Parameter :: ppx_atm_cuall=18      ! All C-grid (u)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_cvall=19      ! All C-grid (v)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_compressed=21 ! Compressed land
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_ozone=22      ! Field on ozone
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_river=23      ! River routing
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_rim=25        ! Rim type field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_theta=26  ! All T points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_u=27      ! All u points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_v=28      ! All v points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_orog=29   ! Orography field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_ocn_nonstd=30     ! Non-standard ocean
                                                  ! grid
      Integer, Parameter :: ppx_ocn_tcomp=31      ! Compressed T points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_ucomp=32      ! Compressed u points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_tall=36       ! All T points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_uall=37       ! All u points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_cuall=38      ! All C-grid (u)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_cvall=39      ! All C-grid (v)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_tfield=41     ! All non-cyclic T
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_ufield=42     ! All non-cyclic u
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_tzonal=43     ! Zonal n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_uzonal=44     ! Zonal n-c field at
                                                  ! u points (ocean)
      Integer, Parameter :: ppx_ocn_tmerid=45     ! Merid n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_umerid=46     ! Merid n-c field at
                                                  ! u points  (ocean)
      Integer, Parameter :: ppx_ocn_scalar=47     ! Scalar (ocean)
      Integer, Parameter :: ppx_ocn_rim=51        ! Rim type field
                                                  ! (LAM BCs ocean)
      Integer, Parameter :: ppx_ocn_lbc_theta=52  ! Ocean rim fields
      Integer, Parameter :: ppx_ocn_lbc_u=53      ! on T & U grids
      Integer, Parameter :: ppx_wam_all=60        ! All points (wave
                                                  ! model)
      Integer, Parameter :: ppx_wam_sea=62        ! Sea points only
                                                  ! (wave model)
      Integer, Parameter :: ppx_wam_rim=65        ! Rim type field
                                                  ! (LAM BCs wave)

!--------------------------------------------------------------------
! Valid rotation type codes
!--------------------------------------------------------------------
      Integer, Parameter :: ppx_unrotated=0       ! Unrotated output
                                                  ! field
      Integer, Parameter :: ppx_elf_rotated=1     ! Rotated ELF field

!-------------------------------------------------------------------
! Valid level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_full_level=1      ! Model full level
      Integer, Parameter :: ppx_half_level=2      ! Model half level
      Integer, Parameter :: ppx_rho_level=1       ! Model rho level
      Integer, Parameter :: ppx_theta_level=2     ! Model theta level
      Integer, Parameter :: ppx_single_level=5    ! Model single level
      Integer, Parameter :: ppx_soil_level=6      ! Deep Soil level

!-------------------------------------------------------------------
! Valid data type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_type_real=1       ! Real data type
      Integer, Parameter :: ppx_type_int=2        ! Integer data type
      Integer, Parameter :: ppx_type_log=3        ! Logical data type

!-------------------------------------------------------------------
! Valid meto8 level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_meto8_surf=9999   ! MetO8 surface type
                                                  ! code

!-------------------------------------------------------------------
! Valid dump packing codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_pack_off=0        ! Field not packed
                                                  ! (ie. 64 bit)
      Integer, Parameter :: ppx_pack_32=-1        ! Field packed to
                                                  ! 32 bit in  dump
      Integer, Parameter :: ppx_pack_wgdos=1      ! Field packed by
                                                  ! WGDOS method
      Integer, Parameter :: ppx_pack_cfi1=11      ! Field packed using
                                                  ! CFI1  (ocean)

!-------------------------------------------------------------------
! Add valid lbvc codes referenced in model (pp header output labels)
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_lbvc_height  =  1 ! height
      Integer, Parameter :: ppx_lbvc_depth   =  2 ! depth (ocean)
      Integer, Parameter :: ppx_lbvc_pressure=  8 ! pressure
      Integer, Parameter :: ppx_lbvc_theta   = 19 ! potential T
      Integer, Parameter :: ppx_lbvc_hybrid  = 65 ! hybrid height(atmos)
      Integer, Parameter :: ppx_lbvc_PV      = 82 ! potential vorticity
      Integer, Parameter :: ppx_lbvc_surface =129 ! surface
! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      INTEGER :: MOS_OUTPUT_LENGTH
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP,MOS_OUTPUT_LENGTH

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
      INTEGER:: MOS_MASK(MOS_MASK_LEN)
! TYPSTS end
!------------------------ nstypes.h ----------------------------------
!jhan:further renovation of ths file may be necessary params are dependent on dataset
!jhan: ALSO nstypes_cable.h should be unecessary nsoil/soil is only difference
      !--- Number of non-vegetation surface types
      Integer, Parameter :: NNVG  = 4

      !--- Number of plant functional types.
      Integer, Parameter :: NPFT  = 13
      
      !--- Number of surface types.
      Integer, Parameter :: NTYPE =17 
      
      !--- Index of the surface type 'Soil'
      !Integer, Parameter :: SOIL  = 16 
      !dhb599, 20110615: change made as per Peter Vohralik, item 1:
      Integer, Parameter :: SOIL  = 14

!--- Land surface types :
!--- original veg. tiles 
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!--- for testing these tiles are set = 1:5 
!     6 - Broadleaf Tree
!     7 - Needleleaf Tree
!     8 - C3 Grass
!     9 - C4 Grass
!    10 - Shrub
!--- for testing these tiles are set = 0
!    11 - 0 
!    11 - 0
!    11 - 0
!--- original non-veg tiles moved to these indices
!     14 - Urban
!     15 - Water
!     16 - Soil
!     17 - Ice


! River routing
! TYPATCPL start
! Description: Gridline coordinates for interpolation and area-averaging
! between atmosphere and river-routing grids (Part of TYPAOCPL.hdk)
! Author: C.Bunton 28.02.03
!
! History:
! Version  Date    Comment
!  5.5  28/02/03   Original code. C. Bunton
!
!
      REAL :: XPA(AOCPL_ROW_LENGTH+1)  ! Atmosphere TP longitude coordina
      REAL :: XUA(0:AOCPL_ROW_LENGTH)  ! Atmosphere U longitude coordinat
      REAL :: XVA(AOCPL_ROW_LENGTH+1)  ! Atmosphere V longitude coordinat
      REAL :: YPA(AOCPL_P_ROWS)        ! Atmosphere TP latitude coordinat
      REAL :: YUA(AOCPL_P_ROWS)        ! Atmosphere U latitude coordinate
      REAL :: YVA(0:AOCPL_P_ROWS)      ! Atmosphere V latitude coordinate
! TYPATCPL end
      INTEGER                                                           &
     &  G_P_FIELD                                                       &
                              ! IN : global horiz domain for atmos
     &, G_R_FIELD             ! IN : global horiz domain for rivers
      INTEGER :: obs_flag_len,obs_len
      INTEGER :: OBS_FLAG(obs_flag_len)
      REAL    :: OBS(obs_len)
!
! CHSUNITS define the number of i/o units
!
!  Author : R A Stratton
!
!  Model            Modification history:
! version  date
!   3.1  03/02/93   Introduced at version 3.1
!   4.1  21/02/96   Increase no.of i/o units to accommodate wave
!                   sub-model.  RTHBarnes.
!   5.2  21/08/00   Add an extra op macro for VAR plus 1 user pp
!                   output stream. R Rawlins
!   6.2  19/01/06   Increased NUNITS to 152 to accomodate extra
!                   diagnostic
!
! Project task:
!
!  Documentation:  Unified Model Documentation Paper
!                  H- History Bricks
!
! ---------------------------------------------------------------

      ! These values must be consistent with OUTFILE_S, OUTFILE_L
      ! and OUTFILE_E in file VERSION.
      INTEGER,PARAMETER::NUNITS=161   ! No. of I/O units
      ! length of most unit no arrays
      INTEGER,PARAMETER::NUNITS_LEN=NUNITS-19

      ! The above parameter statements must not be altered without
      ! considering the effect on the following HISTORY files CHISTO,
      ! CLFHIST and IHISTO.
      ! This file must always preceed the above history file
      ! New file environment variable names may need to be added to
      ! CLFHIST and/or CENVIRDT (usually both) depending on manner of
      ! I/O.
! CHSUNITS end
!LL  Comdeck: CCONTROL -------------------------------------------------
!LL
!LL  Purpose: COMMON block for top level switches and 2nd level switches
!LL           needed by the top level (C0) and 2nd level routines, but
!LL           not held in the history COMMON block.
!LL
!LL Logical components covered :
!LL
!LL External documentation: Unified Model documentation paper No
!LL                         Version
!LL
!LLEND ---------------------------------------------------------------

! ----------------------- Comdeck: CNTLALL  ----------------------------
! Description: COMDECK defining Control variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.0  25/10/95  Add user switch CONTROL_RESUBMIT. RTHBarnes
!  4.4  28/07/97  Add user switch LCLIMREALYR. M Gallani
!  4.4  11/10/97  Add logical switch L_AO_D1_MEMORY. D. Robinson.
!  5.2  14/11/00  Enable Ocean Run Length Encoding. Ian Edmond
!  5.3  25/09/01  Add switch L_IO_Timer. P.Selwood.
!  5.3  18/09/01  Add FT_LASTSTEP. David Baker
!  5.4  17/09/02  Num_ALBCs and ALBC2_StartTime_mins added. Adam Clayton
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
      ! Array holding original data time (prior to assimilation)
      INTEGER:: MODEL_BASIS_TIME(6)

      ! Model analysis time in hours since Basis Time
      ! UM6.5 - Replace MODEL_ANALYSIS_HRS by MODEL_ANALYSIS_MINS 
      !         MODEL_ANALYSIS_HRS changed to REAL
      REAL :: MODEL_ANALYSIS_HRS
      INTEGER :: MODEL_ANALYSIS_MINS

      INTEGER:: MODEL_HRS_PER_GROUP! No. of hours in coupling period
      INTEGER:: NCPU               ! No of CPUs assigned to the program
      INTEGER:: ANCIL_REFTIME(6)   ! Ref. time for updating ancillaries
      INTEGER:: FT_PLOTSEL(60:69)  ! interval for plotting pp file
      INTEGER:: RUN_TARGET_END(6)  ! Target end time for this run

      INTEGER:: Num_ALBCs            ! Number of atmos boundary files
      INTEGER:: ALBC2_StartTime_mins ! VT of first block of data in 2nd
                                     ! atmos boundary file, in minutes
                                     ! from start of run

      ! Increment to be added on each resubmission of the job.
      INTEGER:: RUN_RESUBMIT_INC(6)

      ! Number of field headers reserved for non-mean PPfiles on each
      ! unit
      INTEGER:: PP_LEN2_LOOK(20:NUNITS)

      ! Internally defined PP packing code
      INTEGER:: PP_PACK_CODE(20:NUNITS)

      ! Frequency of initialisation of FTunit
      INTEGER:: FT_STEPS(20:NUNITS)
      INTEGER :: FT_FIRSTSTEP(20:NUNITS)  ! ... starting at step number .
      INTEGER :: FT_LASTSTEP(20:NUNITS)    ! ... ending at step number ..
      LOGICAL:: LATMOSNEXT,LOCEANNEXT ! Flags to select atmosphere/ocean
      LOGICAL:: LPP                   ! Activate PPCTL
      LOGICAL:: LPP_SELECT(20:NUNITS) ! Activate PP init on unit
      LOGICAL:: LDUMP                 ! Activate DUMPCTL
      LOGICAL:: LMEAN                 ! Activate MEANCTL
      LOGICAL:: LHISTORY              ! Update TEMP history file
      LOGICAL:: LPRINT                ! Activate PRINTCTL
      LOGICAL:: LINTERFACE            ! Activate GEN_INTF
      LOGICAL:: LEXIT                 ! Activate EXITCHEK
      LOGICAL:: LJOBRELEASE           ! Activate JOBCTL
      LOGICAL:: LMEANPR(4)            ! Select printed diags from means
      LOGICAL:: LANCILLARY            ! Activate UP_ANCIL
      LOGICAL:: LBOUNDARY             ! Activate UP_BOUND
      LOGICAL:: LASSIMILATION         ! Activate assimilation
      LOGICAL:: LCAL360               ! 360-day calendar
      LOGICAL:: LTIMER                ! Activate detailed TIMER routine
      LOGICAL:: L_AO_D1_MEMORY  ! T : D1 copied to memory for AO coupling
      LOGICAL:: LCLIMREALYR           ! Real-period climate means
      LOGICAL:: LRLE                  ! Indicates Run Length Encoding
      LOGICAL :: L_IO_TIMER              ! Activate IO Timer.

      CHARACTER*4  EXPT_ID          ! Unique alphanumeric serial number
!                                   ! associated with model
!                                   ! (Non-Operational expts)
!                                   !
!                                   ! Operational run name
!                                   ! (Operational expts)
      CHARACTER*8  EXPT_ALIAS       ! Non unique user defined expt name
      CHARACTER*1  JOB_ID           ! Unique alphanumeric job identifier
!                                   ! used for networking
      CHARACTER*4  EXPT_ID_IN       ! Experiment ID of driving model if
!                                   ! limited-area run
      CHARACTER(LEN=1) :: JOB_ID_IN        ! Job ID of driving model if
!                                   ! limited-area run
      CHARACTER*14 MODEL_STATUS     ! Operational or NonOperational
      CHARACTER*14 MODEL_ASSIM_MODE ! Atmosphere,Ocean,Coupled or None
      CHARACTER*17 TIME_CONVENTION  ! Relative, Timestep, Absolute_long,
!                                    Absolute_standard or Absolute_short
      CHARACTER*1  FT_WSSEND(60:69) ! "Y" if file to be sent to HP
!
      CHARACTER*1 TYPE_LETTER_1(20:NUNITS) ! File type letter #1
      CHARACTER*1 TYPE_LETTER_2(20:NUNITS) ! File type letter #2
      CHARACTER*1 TYPE_LETTER_3(20:NUNITS) ! File type letter #3
!
      CHARACTER*1  FT_INPUT (20:NUNITS) ! "Y" if input file on unit
      CHARACTER*1  FT_OUTPUT(20:NUNITS) ! "Y" if output file on unit
      CHARACTER*1  FT_SELECT(20:NUNITS) ! "Y" if file selected for post
!                                          processing request.
      CHARACTER*1  FT_ARCHSEL(20:NUNITS) ! "Y" if file to be archived.
!
      CHARACTER*10 RUN_ASSIM_MODE      ! cf MODEL_ASSIM_MODE (Oper use)
      CHARACTER*1  CONTROL_RESUBMIT    ! User flag for auto resubmit

      NAMELIST / NLSTCALL /                                             &
     &  MODEL_BASIS_TIME, MODEL_ANALYSIS_MINS,                          &
     &  MODEL_HRS_PER_GROUP,                                            &
     &  NCPU, ANCIL_REFTIME, FT_PLOTSEL, RUN_TARGET_END,                &
     &  Num_ALBCs, ALBC2_StartTime_mins,                                &
     &  RUN_RESUBMIT_INC, PP_LEN2_LOOK, PP_PACK_CODE,                   &
     &  FT_STEPS, FT_FIRSTSTEP, FT_LASTSTEP,                            &
     &  LATMOSNEXT, LOCEANNEXT, LPP, LPP_SELECT, LDUMP, LMEAN,          &
     &  LHISTORY, LPRINT, LINTERFACE, LEXIT, LJOBRELEASE,               &
     &  LMEANPR, LANCILLARY, LBOUNDARY, LASSIMILATION,                  &
     &  LCAL360, LTIMER, L_AO_D1_MEMORY, LRLE,                          &
     &  LCLIMREALYR, L_IO_TIMER,                                        &
     &  EXPT_ID, JOB_ID, EXPT_ID_IN, JOB_ID_IN,                         &
     &  EXPT_ALIAS, MODEL_STATUS, MODEL_ASSIM_MODE,                     &
     &  TIME_CONVENTION, FT_WSSEND,                                     &
     &  TYPE_LETTER_1, TYPE_LETTER_2, TYPE_LETTER_3,                    &
     &  FT_INPUT, FT_OUTPUT, FT_SELECT, FT_ARCHSEL,                     &
     &  RUN_ASSIM_MODE, CONTROL_RESUBMIT
      COMMON / CNTLCALL /                                               &
     &  MODEL_BASIS_TIME, MODEL_ANALYSIS_MINS,                          &
     &  MODEL_HRS_PER_GROUP,                                            &
     &  NCPU, ANCIL_REFTIME, FT_PLOTSEL, RUN_TARGET_END,                &
     &  Num_ALBCs, ALBC2_StartTime_mins,                                &
     &  RUN_RESUBMIT_INC, PP_LEN2_LOOK, PP_PACK_CODE,                   &
     &  FT_STEPS, FT_FIRSTSTEP, FT_LASTSTEP,                            &
     &  LATMOSNEXT, LOCEANNEXT, LPP, LPP_SELECT, LDUMP, LMEAN,          &
     &  LHISTORY, LPRINT, LINTERFACE, LEXIT, LJOBRELEASE,               &
     &  LMEANPR, LANCILLARY, LBOUNDARY, LASSIMILATION,                  &
     &  LCAL360, LTIMER, L_AO_D1_MEMORY, LRLE,                          &
     &  LCLIMREALYR, L_IO_TIMER,                                        &
! Character variables at the end of the common block
     &  EXPT_ID, JOB_ID, EXPT_ID_IN, JOB_ID_IN,                         &
     &  EXPT_ALIAS, MODEL_STATUS, MODEL_ASSIM_MODE,                     &
     &  TIME_CONVENTION, FT_WSSEND,                                     &
     &  TYPE_LETTER_1, TYPE_LETTER_2, TYPE_LETTER_3,                    &
     &  FT_INPUT, FT_OUTPUT, FT_SELECT, FT_ARCHSEL,                     &
     &  RUN_ASSIM_MODE, CONTROL_RESUBMIT
! ----------------------- Comdeck: CNTLGEN  ----------------------------
! Description: COMDECK defining Control variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  28/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.0   3/11/95  Move character array MEANSim to the end of the
!                 common block to ensure that it starts correctly on a
!                 word boundary. [No problem is apparent on the Cray
!                 if N_INTERNAL_MODEL_MAX is an even no.]
!                 Rick Rawlins
!  4.1  03/04/96  Add new array DUMP_PACKim. D. Robinson
!  4.5  10/11/98  Increase number of dumps allowed at irregular
!                 timesteps from 10 to 40: Move lengths into
!                 CNTLGEN. R Rawlins
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!

      ! Max no. of irregular times for dumps
      INTEGER, PARAMETER :: DUMPTIMES_LEN1 = 40

      ! No. of areas of zonal mean prints
      INTEGER, PARAMETER :: PRINTFREQ_LEN1 = 5

      ! No. of time intervals for climate meaning
      INTEGER, PARAMETER :: MEANFREQ_LEN1 = 4

      ! Max no. of irregular times for job release
      INTEGER, PARAMETER :: JOBREL_LEN1 = 10

      INTEGER:: STEPS_PER_PERIODim(N_INTERNAL_MODEL_MAX)
      INTEGER:: SECS_PER_PERIODim(N_INTERNAL_MODEL_MAX)

      ! Number of advection timesteps between checks for model exit
      INTEGER:: EXITFREQim(N_INTERNAL_MODEL_MAX)

      ! Number of steps between atmosphere restart dumps
      INTEGER:: DUMPFREQim(N_INTERNAL_MODEL_MAX)

      ! Archiving frequency  for atmos dumps
      INTEGER:: ARCHDUMP_FREQim(N_INTERNAL_MODEL_MAX)

      ! Timesteps (from start of run) at which restart dumps are written
      INTEGER:: DUMPTIMESim(DUMPTIMES_LEN1,N_INTERNAL_MODEL_MAX)

      ! Indicators for mean dump frequency
      INTEGER:: MEANFREQim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Switches for mean dump arch.
      INTEGER:: MEANARCHim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! PP field selectors
      INTEGER:: PPSELECTim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Switches for pp field archive
      INTEGER:: ARCHPPSELim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Switches for chart plotting
      INTEGER:: PLOTSELim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Number of field headers to reserve for internal model mean
      ! PPfiles
      INTEGER:: PP_LEN2_MEANim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Reference time for production of means
      INTEGER:: MEAN_REFTIMEim(6,N_INTERNAL_MODEL_MAX)

      ! Indicators of zonal mean print frequency
      INTEGER:: PRINTFREQim(PRINTFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Step numbers  at which to release user-specified scripts

      INTEGER:: JOBREL_STEPim(JOBREL_LEN1,N_INTERNAL_MODEL_MAX)

      ! Offset for dump archiving
      INTEGER:: ARCHDUMP_OFFSETim(N_INTERNAL_MODEL_MAX)

      ! Unit reserved for mean PPs
      INTEGER:: FT_MEANim(N_INTERNAL_MODEL_MAX)

      ! Packing indicator for dumps
      INTEGER:: DUMP_PACKim(N_INTERNAL_MODEL_MAX)

      ! "Y" if mean file to be sent to HP
      CHARACTER(LEN=1) :: MEANWSim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      LOGICAL:: LLBOUTim(N_INTERNAL_MODEL_MAX)  ! Lateral b.c.'s
      LOGICAL:: LANCILim(N_INTERNAL_MODEL_MAX)  ! Ancillary files

      NAMELIST / NLSTCGEN /                                             &
     &  STEPS_PER_PERIODim, SECS_PER_PERIODim,                          &
     &  EXITFREQim, DUMPFREQim,                                         &
     &  ARCHDUMP_FREQim, DUMPTIMESim, PPSELECTim, PLOTSELim,            &
     &  ARCHPPSELim, MEANARCHim, MEANFREQim, MEAN_REFTIMEim,            &
     &  PRINTFREQim,  JOBREL_STEPim, ARCHDUMP_OFFSETim, PP_LEN2_MEANim, &
     &  FT_MEANim,                                                      &
     &  DUMP_PACKim,                                                    &
     & MEANWSim, LLBOUTim, LANCILim

      COMMON / CNTLCGEN /                                               &
     &  STEPS_PER_PERIODim, SECS_PER_PERIODim,                          &
     &  EXITFREQim, DUMPFREQim,                                         &
     &  ARCHDUMP_FREQim, DUMPTIMESim, PPSELECTim, PLOTSELim,            &
     &  ARCHPPSELim, MEANARCHim, MEANFREQim, MEAN_REFTIMEim,            &
     &  PRINTFREQim,  JOBREL_STEPim, ARCHDUMP_OFFSETim, PP_LEN2_MEANim, &
     &  FT_MEANim,                                                      &
     &  DUMP_PACKim,                                                    &
     &  LLBOUTim, LANCILim,                                             &
     &  MEANWSim

! CNTLGEN end
! ----------------------- header file: CNTLATM  -----------------------
! Description: Control variables for the Atmosphere model (read only).
!              Contains logical switches for science options needed
!              for addressing calculations and intermediate control.
!              Predominantly used for holding logical flags set up by
!              the UMUI - read in by a namelist, but can also hold
!              derived control flags set by the model.
!              [Note that CRUNTIMC holds accompanying run-time
!              constants.]
!
! Author : R. Rawlins
!
! History:
! Version  Date      Comment.
!  5.0 20/04/99  Extensive revision to earlier comdeck for conversion
!                to C-P 'C' dynamics grid. R. Rawlins
!  5.1  4/02/00  Restore energy correction switches removed at UM5.0
!                Also added additional switches for mass and moisture
!                R A Stratton.
!  5.1 13/04/00  IAU control moved to CTFilt. Adam Clayton
!  5.2 22/08/00  Reinstate Murk and Tracer switches. P.Selwood.
!  5.2 15/11/00  Reintroduce logicals for MOSES 2 and triffid. M.Best.
!  5.3   12/10/01   Remove hard-wired L_trivial_trigs.
!                   put problem_number in CNTLATM.    Terry Davies
!  5.3 27/07/01  Add logical switch L_MURK_RAD   S. Cusack
!  5.3    06/01  Introduce and hardwire logical for setting of
!                leads temperature.                         Nic Gedney
!  5.3 15/10/01  Added L_USE_METHOX switch. David Jackson
!  5.3 29/08/01  Sulphate and sea-salt logicals split into separate
!                versions for different processes.  A. Jones
!  5.3 09/04/01  Add logical switch for spatial degradation of
!                E-S radiation calculations.             S. Cusack
!
!  5.3 19/06/01   Stuff to handle tropopause-based ozone added
!                 -- see Radiation block              Dave Tan
!  5.4 15/08/02   Reconfiguration use changes. P.Selwood.
!  5.4  2/08/02  Add logical switch for PC2 cloud and condensation
!                scheme.                              Damian Wilson
!  5.4 24/10/02  Moved L_GWD and L_USE_USSP from CRUNTIMC. S. Webster
!  5.4 11/03/02  Remove comment lines from after #include
!                                                 S. Carroll
!  5.5 06/02/03  River routing support. P.Selwood.
!  5.5 05/02/03  Add logicals for biomass aerosol scheme     P Davison

!  5.5 17/02/03  Add two large-scale hydrology logicals.
!                L_TOP and L_PDM.                  Nic Gedney
!  5.5 08/01/03  Remove L_3DVAR, L_4DVAR, L_3DVAR_BG, LSINGLE_HYDROL,
!                L_H2_SULPH and L_LSPICE_BDY. D Robinson
!  5.5 21/01/03  Move L_USE_TPPS_OZONE to NLSTCATM namelist. D Robinson
!  5.5 13/03/03  Move I_OZONE_INT here from CRUNTIMC.
!                                                  Jean-Claude Thelen
!  5.5 03/02/03  Include L_mcr logicals.    R.M.Forbes
!  5.5 19/02/03  Remove redundant L_BL_LSPICE and L_RMBL   A.Lock
!  6.0 30/07/03  Include l_pc2_lbc for cloud frac. lbcs. Damian Wilson
!  6.0 19/08/03  Add runtime controls for 4 physics sections
!  6.1  02/08/04  Add logicals for stochem coupling. R Barnes
!  6.1 13/05/04  Add super_array_size variable                 A.Malcolm
!  6.1 07/04/04  Add logical for autoconversion de-biasing.   A. Jones
!  6.1 07/04/04  Add logicals for interactive DMS emissions.  A. Jones
!  6.2 25/01/06  Add iteration count logical test variable    T.Edwards
!  6.2 15/11/05  Add logical for aerosol optical depth     N. Bellouin
!  6.2 23/11/05  Add logical for RH and hygroscopic aerosols N.Bellouin
!  6.2 01/03/06  Add L_UPPER_STRATQ switch - David Jackson
!  6.2 06/01/06  Add logicals for seaice albedo options. J.Ridley
!  6.2 23/02/06  Add logicals for UKCA sub-model.  F.O'Connor
!  6.1 07/04/04  Add logical for RothC temperature function.  C.D. Jones
!  6.2 21/07/05  Add moisture_array_size variable              A.Malcolm
!  6.2 01/10/05  Include L_murk_lbc logical.    R.M.Forbes
!  6.2 07/11/05  Add L_MOD_BARKER_ALBEDO and L_USE_SPEC_SEA
!                to NLSTCATM.        James Manners
!  6.2 25/01/06  Add L_INHOM_CLOUD to NLSTCATM.    James Manners
!  6.2 09/03/06  Add logicals for inland basin rerouting. P.Falloon
!  6.2 24/02/06  Add logical for 10m windspeed pass A2O  J.Gunson

!  6.2 24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!  6.2 07/11/05 Add L_use_orog_corr to NLSTCATM.     James Manners
!  6.2 24/02/06  Add logical to allow DMS flux from ocean model J.Gunson
!  6.4 19/01/07 Removed a comment relating to A05_3c scheme. R A Stratton
!  6.4 08/01/06 Include Brooks cloud fraction logicals. Damian Wilson
!  6.4 10/01/07 Include mixing ratio control logical. Damian Wilson
!  6.6.2 10/06/09 Logicals for reading solar/volcanic forcing. Chris Jones
!------------------   General:  -------------------------------------
      INTEGER Model_domain        ! Domain of atmosphere model:
!                                   global,LAM,cyclic LAM,single column
! Model_domain meaningful names
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end

      INTEGER :: problem_number      ! type of problem to be solved

      INTEGER MAXSECTS            ! Max. no. of code sections
      PARAMETER (MAXSECTS=99)
      CHARACTER*3 H_SECT(0:MAXSECTS) ! Array of code section-versions

      ! Physics:   ------------------------------------
      LOGICAL :: l_ssice_albedo     ! Sea-ice albedo affected by snow

      LOGICAL :: l_sice_heatflux    ! Semi-impl update of sea-ice temp

      LOGICAL :: l_sice_meltponds ! Sea-ice albedo affected by meltponds

      LOGICAL :: l_sice_scattering  ! Sea-ice albedo affected scattering

      LOGICAL :: l_sice_hadgem1a ! HadGEM1 sea-ice albedo bug corrected

      LOGICAL :: L_NEG_TSTAR        ! Test for -ve surface temperature.
!
      ! Use sulphate aerosol in microphysics
      LOGICAL :: L_USE_SULPHATE_AUTOCONV

      ! Use sea-salt aerosol in microphysics
      LOGICAL :: L_USE_SEASALT_AUTOCONV

      ! Use soot aerosol in microphysics
      LOGICAL :: L_USE_SOOT_AUTOCONV

      ! Use biomass aerosol in microphysics
      LOGICAL :: L_USE_BMASS_AUTOCONV
      
      ! Use fossil-fuel organic carbon in microphysics
      LOGICAL :: L_USE_OCFF_AUTOCONV
      
      ! Use autoconversion de-biasing scheme in microphysics
      LOGICAL :: L_AUTO_DEBIAS
      ! Use sulphate aerosol no. in S-cycle
      LOGICAL :: L_USE_SULPHATE_SULPC

      ! Use sea-salt aerosol no. in S-cycle
      LOGICAL :: L_USE_SEASALT_SULPC

      ! Use soot aerosol no. in S-cycle
      LOGICAL :: L_USE_SOOT_SULPC

      ! Use biomass aerosol no. in S-cycle
      LOGICAL :: L_USE_BMASS_SULPC
      
      ! Use fossil-organic carbon aerosol no. in S-cycle
      LOGICAL :: L_USE_OCFF_SULPC
      
      ! Energy correction:
      LOGICAL :: L_emcorr    ! T: turns on energy correction code
      LOGICAL :: LMASS_corr  ! T: Apply mass correction
      LOGICAL :: LQT_corr    ! T: Apply total moisture correction
      LOGICAL :: LEMQ_print  ! T: Print additional info from em code
      LOGICAL :: LENERGY     ! T: if timestep to cal energy correction
      LOGICAL :: LFLUX_RESET ! T: if timestep to reset flux array in D1

      ! number of model timesteps per energy correction period.
      INTEGER :: A_ENERGYSTEPS

      ! Radiation:

      LOGICAL :: L_radiation     !  F: Turns off radiation code
      LOGICAL :: L_MICROPHY           !  Microphysics in sw rad scheme

!     Use mineral dust in radiation calculations
      LOGICAL :: L_USE_DUST

!     Use biogenic aerosol in radiation code
      LOGICAL :: L_USE_BIOGENIC

      ! Use SO4 aerosol from sulphur cycle for direct/indirect effect
      ! in radiation, the latter for both SW and LW.
      LOGICAL :: L_USE_SULPC_DIRECT
      LOGICAL :: L_USE_SULPC_INDIRECT_SW
      LOGICAL :: L_USE_SULPC_INDIRECT_LW
      ! Indirect radiative effect of sea-salt
      LOGICAL :: L_USE_SEASALT_INDIRECT

      ! Direct radiative effect of sea-salt
      LOGICAL :: L_USE_SEASALT_DIRECT
      LOGICAL :: L_USE_SOOT_DIRECT  ! direct radiative effects of soot
      LOGICAL :: L_USE_SOOT_INDIRECT  ! indirect effects of soot
      ! Use biomass aerosol for direct/indirect effect in radiation.
      LOGICAL :: L_USE_BMASS_DIRECT
      LOGICAL :: L_USE_BMASS_INDIRECT
      
      ! Use fossil-fuel organic carbon aerosol for direct/indirect
      ! effect in radiation
      LOGICAL :: L_USE_OCFF_DIRECT
      LOGICAL :: L_USE_OCFF_INDIRECT
      
!     Use aerosol climatologies in radiation instead of prognostic variables
!     Set on a species by species basis
      LOGICAL :: L_USE_ARCLBIOM   ! biomass burning aerosol
      LOGICAL :: L_USE_ARCLBLCK   ! black carbon
      LOGICAL :: L_USE_ARCLSSLT   ! sea salt
      LOGICAL :: L_USE_ARCLSULP   ! sulpahtes
      LOGICAL :: L_USE_ARCLDUST   ! mineral dust
      LOGICAL :: L_USE_ARCLOCFF   ! organic carbon (fossil fuel)
      LOGICAL :: L_USE_ARCLDLTA   ! delta aerosol

      ! Aerosol optical depth diagnostic was requested
      LOGICAL :: L_USE_AOD

      ! Clear-sky relative humidity is to be used instead of
      ! grid-box mean RH, for hygroscopic aerosols
      LOGICAL :: L_USE_CLEARRH

      ! controls the use of spatial degradation of radiation calc.
      LOGICAL :: L_rad_deg
      LOGICAL :: L_CTILE       ! Switch for coastal tiling.
!                              ! If TRUE then land and sea can
!                              ! coexist in the same gridbox.
!                              ! If FALSE the land fraction
!                              ! must be equal to 0 to 1
      LOGICAL :: L_TOP         ! If TRUE then TOPMODEL-based
!                              ! hydrology scheme.
      LOGICAL :: L_PDM         ! If TRUE then PDM-based
!                              ! hydrology scheme.
      LOGICAL :: L_SOIL_SAT_DOWN ! If TRUE then super-saturated soil
!                                ! moisture moves downward, else upward

      INTEGER :: H_SWBANDS    ! Number of shortwave radiation bands
      INTEGER :: H_LWBANDS    ! Number of longwave radiation bands
      INTEGER :: A_SW_RADSTEP ! Number of advection steps per SW step
      INTEGER :: A_LW_RADSTEP ! Number of advection steps per LW step
! Number of advection steps per prognostic/diagnostic SW and LW step.
!'Prognostic' and 'Diagnostic' refer to the frequency of the calls
! to radiation code in the Unified Model.
! In the case of time stepping prognostic and diagnostic refer to the
! slow and fast radiative timestep respectively. In the case of radiative
! forcing they refer to the prognostic and diagnostic calls to radiation.
      INTEGER :: A_SW_RADSTEP_DIAG
! Number of advection steps per 'fast' SW step (3C)
      INTEGER :: A_LW_RADSTEP_DIAG
! Number of advection steps per 'fast' LW step (3C)
      INTEGER :: A_SW_RADSTEP_PROG
! Number of advection steps per 'slow' LW step (3C)
      INTEGER :: A_LW_RADSTEP_PROG
! Number of advection steps per 'slow' LW step (3C)

      INTEGER :: i_ozone_int  ! Option for interpolation of ozone

      ! Cloud:

      LOGICAL:: L_CLD_AREA           ! controls cloud area parametriz.
      LOGICAL:: L_ACF_CUSACK         ! ... whether to have Cusack
      LOGICAL:: L_ACF_BROOKS         ! ... or Brooks
      LOGICAL:: L_PC2                ! controls PC2 cloud scheme
      LOGICAL:: L_PC2_RESET          ! run PC2 scheme diagnostically
      LOGICAL:: L_PC2_LBC            ! LBC's contain cloud fractions
      LOGICAL:: L_PC2_DIAG_SH        ! Use diagnostic convective shallow cloud
                                     ! in PC2.

      ! Assimilation:

      ! Switches for assm mode
      LOGICAL:: L_AC

      ! T: Use V_INT_TP to output Temp on model levels.
      LOGICAL:: L_VINT_TP

      ! UM6.5 - MODEL_ANALYSIS_HRS replaced by MODEL_ANALYSIS_MINS - 
      !         change A_ASSIM_START_HR and A_ASSIM_END_HR to 
      !         A_ASSIM_START_MIN, A_ASSIM_END_MIN
      !         so that all three variables  have the same flexibility
      ! Time at which data assimilation starts (Hours after Basis Time)
      INTEGER :: A_ASSIM_START_MIN
      ! Time at which data assimilation  ends
      INTEGER :: A_ASSIM_END_MIN

      ! Switch for assimilation mode
      CHARACTER(LEN=5) :: A_ASSIM_MODE
      
      !---  PMSL diagnostic  ---
      ! Orographic height threshold for new pmsl calculation
      ! from ATMOS_STASH_Misc in UMUI for Calc_NPMSL routine
      REAL :: NPMSL_HEIGHT

      ! Switch for interpolated winds in lbcs
      LOGICAL :: L_int_uvw_lbc  ! .true. for advecting winds
                                ! interpolated in boundary zone
      
      !---  Tracers ---
      ! Aerosol

      LOGICAL:: L_MURK          !      :Total aerosol field
      LOGICAL:: L_MURK_ADVECT   !      :Aerosol advection
      LOGICAL:: L_MURK_SOURCE   !Bndry :Aerosol source & sink terms
      LOGICAL:: L_MURK_BDRY     !Layer :UK Mes bndry model
      LOGICAL:: L_BL_TRACER_MIX !model :Bndry layer tracer mixing
      LOGICAL :: L_MURK_RAD
      LOGICAL :: L_murk_lbc    !  Murk aerosol lbcs active

      ! For Aero_Ctl (Sulphur cycle or Soot)

      INTEGER CALL_CHEM_FREQ     !No. times chem called per atm tstep

!     Mineral dust aerosol
      LOGICAL :: L_DUST
!     Use old version of dust_uplift scheme used in CAM NWP models
      LOGICAL :: L_CAM_DUST

      ! Sulphur cycle

      LOGICAL :: L_SULPC_SO2   ! S Cycle: SO2 MMR included
      LOGICAL :: L_SULPC_DMS   ! S Cycle: DMS MMR included
      LOGICAL :: L_SULPC_OZONE ! S Cycle: Ozone included for oxidation 
                               !          of DMS and SO2
      LOGICAL :: L_SULPC_SO2_O3_NONBUFFERED ! S Cycle: SO2+O3 reaction NOT
                                            ! buffered by NH3.
      LOGICAL :: L_SULPC_ONLINE_OXIDANTS ! Sulphur Cycle : Use online
                                         ! oxidants from UKCA
      LOGICAL :: L_SO2_SURFEM  ! SO2 Surface Emissions
      LOGICAL :: L_SO2_HILEM   ! SO2 High Level Emissions
      LOGICAL :: L_SO2_NATEM   ! SO2 Natural Emissions
      LOGICAL :: L_DMS_EM      ! DMS Emissions
      LOGICAL :: L_DMS_EM_INTER      ! Interactive DMS Emissions
      LOGICAL :: L_DMS_Ointer        ! DMS emissions from ocean model
      LOGICAL :: L_DMS_Liss_Merlivat ! Switches to determine which
      LOGICAL :: L_DMS_Wanninkhof    !    scheme to use for interactive
      LOGICAL :: L_DMS_Nightingale   !    sea-air exchange of DMS
      LOGICAL :: L_SULPC_NH3   ! S Cycle: NH3 tracer included
      LOGICAL :: L_NH3_EM      ! S Cycle: NH3 emiss included

      ! Soot cycle

      LOGICAL :: L_SOOT                ! Soot included
      LOGICAL :: L_SOOT_SUREM          ! surface Soot emiss included
      LOGICAL :: L_SOOT_HILEM          ! elevated Soot emiss included

      ! Biomass aerosol

      LOGICAL :: L_BIOMASS             ! Biomass aerosol included
      LOGICAL :: L_BMASS_SUREM         ! Sfc biomass emiss included
      LOGICAL :: L_BMASS_HILEM         ! Elevated bmass emiss included
      
      ! Fossil-fuel organic carbon aerosol
      
      LOGICAL :: L_OCFF                ! OCFF aerosol included
      LOGICAL :: L_OCFF_SUREM          ! Surface OCFF emissions included
      LOGICAL :: L_OCFF_HILEM          ! Elevated OCFF emiss included
      
      ! Carbon cycle

      ! Interactive 3D CO2 field for use with carbon cycle model
      LOGICAL :: L_CO2_INTERACTIVE
      ! Switch for Radiation Interaction with CO2 - kdcorbin, 06/10
      LOGICAL :: L_CO2_RADIATION   
      ! Switch for CABLE - kdcorbin, 03/10
      LOGICAL :: l_cable
      
      ! Switch for calculating CO2/tracer atmospheric mass - kdcorbin, 05/10
      LOGICAL :: L_TRACER_MASS, L_CO2_MASS
      INTEGER :: I_TRACERMASS_START
      ! Switch for running passive tracers using CO2 fluxes - rml, 1/7/13
      LOGICAL :: L_CO2_TRACER

      ! Switch for calculating methane atmospheric loss - kdcorbin, 05/10
      LOGICAL :: L_METHANE_LOSS
      INTEGER :: I_METHANE_TRACERS

      ! Switch for calculating mcf atmospheric/ocean loss - kdcorbin, 05/10
      LOGICAL :: L_MCF_LOSS
      INTEGER :: I_MCF_TRACERNUMBER

      ! Switch for calculating radon decay - kdcorbin, 05/10
      LOGICAL :: L_RADON_DECAY
      INTEGER :: I_RADON_TRACERNUMBER

      ! CO2 Mass - kdcorbin, 05/10
      REAL :: CO2EMITMASS

     ! Tracer Mass - kdcorbin, 05/10
      REAL :: TMASS(21)

      LOGICAL :: L_CO2_EMITS          ! Include surface emissions

      ! 10m windspeed for air/sea gas flux calculations
      LOGICAL :: L_PSSWND          ! pass 10m windspeed to the ocean

      ! Dust deposition for ocean biology
      LOGICAL :: L_DUST2OCN        ! pass dust dep fields to the ocean

      LOGICAL :: L_Q10                  ! control T fn for soil resp

      ! Switch for turning off boundary layer code

      LOGICAL :: L_bl            !  F: Turns off boundary layer code
      ! MOSES II and Triffid logicals--------------------

      LOGICAL :: L_VEG_FRACS          ! Switch for vegetation fractions
      LOGICAL :: L_SNOW_ALBEDO        ! Prognostic snow albedo
      LOGICAL :: L_TRIFFID            ! Switch for interactive veg model
      LOGICAL :: L_PHENOL             ! Switch for leaf phenology

      ! Switch for running TRIFFID in equilibrium mode
      LOGICAL :: L_TRIF_EQ

      ! Switch for starting NRUN mid-way through a TRIFFID calling
      ! period
      LOGICAL :: L_NRUN_MID_TRIF
      LOGICAL :: L_DISTURB      ! Switch for anthropogenic disturbance

      INTEGER :: CAN_MODEL ! Switch for thermal vegetation canopy

      ! Vegetation:

      ! Update frequency for leaf phenology (days)
      INTEGER :: PHENOL_PERIOD

      INTEGER :: TRIFFID_PERIOD ! Update frequency for TRIFFID (days)

      ! Switch for anthropogenic heat source 
      LOGICAL :: l_anthrop_heat_src 

      ! Hardwire until re-assessment of whether these need to be
      ! re-introduced as UMUI-set switches.
! RR old switches needed for addressing but should be defunct?
! RR - can be set with parameters if needed in the interim. ->

      ! Large scale precipitation:

      LOGICAL :: L_rain          !  F: Turns off precipitation code

      ! 'New' cloud/precip microphysics, Defunct, only mixed phase
      ! phys supported
      LOGICAL, PARAMETER :: L_LSPICE    =.true.

      ! Microphysics complexity
      LOGICAL :: L_mcr_qcf2    !  Include second ice variable
      LOGICAL :: L_mcr_qrain   !  Include prognostic rain
      LOGICAL :: L_mcr_qgraup  !  Include prognosic graupel
      LOGICAL :: L_mcr_qcf2_lbc    !  Second ice variable lbcs active
      LOGICAL :: L_mcr_qrain_lbc   !  Prognostic rain lbcs active
      LOGICAL :: L_mcr_qgraup_lbc  !  Prognosic graupel lbcs active

      ! Controls the use of new RHcrit parametrization, option in Sec 9
      ! vn 2A.
      LOGICAL :: L_RHCPT

! Logicals for different radiation packages
      LOGICAL :: L_FORCING
! Calculate radiative forcings (3C)
      LOGICAL :: L_RADIANCE
! Calculate radiances          (3C)
      LOGICAL :: L_TIMESTEP
! Use new timestepping scheme  (3C)
      LOGICAL :: L_WENYI       
! Include Wenyi's pressure & temperature scaling (3A/3C)

      ! Convection:

      LOGICAL :: L_3D_CCA             ! Use 3D conv cloud amount
      LOGICAL :: L_PHASE_LIM          ! Limits phase change of precip
                                      ! in convective downdraught
      LOGICAL :: L_CCRAD              ! Main logical, will remove code
                                      ! connected with CCRad
                                      ! (including bugfixes)
      LOGICAL :: L_3D_CCW             ! Requires l_ccrad=.TRUE.
                                      ! .TRUE. : Radiation to use 3d ccw
                                      ! profile passed to it from
                                      ! convection.
                                      ! .FALSE.: Radiation constructs
                                      ! mean CCW profile from cclwp,ccb
                                      ! and cct as in original.

      ! Timestep frequency for calling convection
      ! Hardwired to calling every timestep
      INTEGER,PARAMETER :: A_CONV_STEP = 1

      ! GWD scheme:
      LOGICAL :: L_GWD        ! Use SSO drag scheme
      LOGICAL :: L_USE_USSP   ! Use spectral GWD scheme

      ! Radiation:

      ! Changes to open sea albedo for HadGEM1
      LOGICAL :: L_MOD_BARKER_ALBEDO ! Modified Barker albedo
      LOGICAL :: L_USE_SPEC_SEA      ! Spectr. dep. sea albedos

      ! Use modulus of fluxes to remove negative effective extinctions
      LOGICAL :: L_MOD_K_FLUX

      ! Fix the selection of fractional sea points in LW radiation
      LOGICAL :: L_CTILE_FIX

      ! Fix instability in quadratic correction to LW source term
      LOGICAL :: L_QUAD_SRC_FIX

      ! Scale the condensed water content to simulate
      ! inhomogeneous clouds
      LOGICAL :: L_INHOM_CLOUD

! Orography correction to SW radiation
      LOGICAL :: L_use_orog_corr    !  Find gradients from mean orog
      LOGICAL :: L_use_grad_corr    !  Use ancillary X & Y gradients

      ! Tropopause-based Ozone Scheme
      LOGICAL :: L_use_tpps_ozone   !  Use TPPS ozone scheme

! Methane oxidation
      REAL    :: Z_TOP
      LOGICAL :: L_USE_METHOX

! STOCHEM coupling to radiation
      LOGICAL :: L_USE_STOCHEM_CH4   ! for methane
      LOGICAL :: L_USE_STOCHEM_O3    ! for ozone
      ! River Routing
      LOGICAL :: L_RIVERS
      LOGICAL :: L_INLAND   ! control rerouting of inland basin water
      REAL    :: RIVER_STEP

      ! Hydrology:

      LOGICAL :: L_hydrology     !  F: Turns off hydrology code

! Max humidity in STRATQ
      LOGICAL :: L_UPPER_STRATQ

      LOGICAL, PARAMETER :: LMOSES        =.TRUE.  ! MOSES hydrology
      LOGICAL :: L_ICOUNT       !  T: Output iteration counts

      ! Mixing ratios:

      Logical :: l_mr_physics1            ! Use mixing ratio in
                                          ! atmos_physics1
      Logical :: l_mr_physics2            ! Use mixing ratio in
                                          ! atmos_physics2
! Stochastic Physics Random Parameters      
      LOGICAL :: L_RPSEED_READ  !  T: Read in previously specified seed
      LOGICAL :: L_RPSEED_WRITE !  T: WRITE out seed


! Ozone tracer as input to radiation scheme      
      LOGICAL :: L_USE_CARIOLLE
      LOGICAL :: L_USE_OZONEINRAD

! RR old switches---------------------------------------- <-

! OASIS coupling
      LOGICAL :: L_OASIS   ! OASIS coupling switch
      LOGICAL :: L_COUPLE_MASTER    ! Couple through master PE
      INTEGER :: OASIS_COUPLE_FREQ  ! Coupling frequency in
                                    ! number of timesteps. 

!     Logicals for UK Chemistry and Aerosols (UKCA) Model

      LOGICAL :: L_ukca           ! True when UKCA is switched on

! Natural climate forcing
      LOGICAL :: L_SCVARY            ! time varying solar forcing
      LOGICAL :: L_VOLCTS            ! time varying volcanic forcing

      COMMON/ CNTLCATM/                                                 &
     &  Model_domain,L_emcorr,                                          &
     &  L_OASIS, OASIS_COUPLE_FREQ, L_COUPLE_MASTER,                    &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                          &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,              &
     &  problem_number,                                                 &
     &  L_SNOW_ALBEDO,l_ssice_albedo,L_MICROPHY,H_SWBANDS,H_LWBANDS,    &
     &  A_SW_RADSTEP,A_LW_RADSTEP,                                      &
     &  A_SW_RADSTEP_DIAG,A_LW_RADSTEP_DIAG,                            &
     &  A_SW_RADSTEP_PROG,A_LW_RADSTEP_PROG,                            &
     &  L_CLD_AREA,L_ACF_CUSACK,L_ACF_BROOKS,L_PC2,L_PC2_RESET,         &
     &  L_PC2_LBC,L_PC2_diag_sh,L_rad_deg, L_AC ,L_VINT_TP,             &
     &  A_ASSIM_START_MIN, A_ASSIM_END_MIN,                             &
     &  NPMSL_HEIGHT,                                                   &
     &  LENERGY, LFLUX_RESET, LMASS_CORR , LQT_CORR, LEMQ_PRINT,        &
     &  A_ENERGYSTEPS,L_GWD,L_USE_USSP,                                 &
     &  L_Murk, L_MURK_ADVECT, L_MURK_SOURCE, L_MURK_BDRY,              &
     &  L_MURK_RAD, L_murk_lbc, L_BL_TRACER_MIX, L_int_uvw_lbc,         &
     &  L_DUST,L_CAM_DUST,L_SULPC_SO2,L_SULPC_DMS,L_SULPC_OZONE,        &
     &  L_SULPC_SO2_O3_NONBUFFERED,L_SO2_SURFEM, L_SO2_HILEM,           &
     &  L_SO2_NATEM,L_DMS_EM, L_DMS_EM_INTER,                           &
     &  L_SULPC_ONLINE_OXIDANTS,                                        &
     &  L_DMS_Ointer,                                                   &
     &  L_DMS_Liss_Merlivat, L_DMS_Wanninkhof, L_DMS_Nightingale,       &
     &  L_SULPC_NH3, L_NH3_EM, L_SOOT, L_SOOT_SUREM, L_SOOT_HILEM,      &
     &  L_BIOMASS, L_BMASS_SUREM, L_BMASS_HILEM,                        &
     &  L_PSSWND, L_DUST2OCN,                                           &
     &  L_RHCPT, L_CCRAD, L_3D_CCA, L_3D_CCW, L_PHASE_LIM,              &
     &  L_CO2_INTERACTIVE,L_CO2_EMITS,l_cable,                          &
! rml 1/7/13
     &  L_CO2_TRACER,                                                   &
!kdcorbin, 08/10
     &  L_CO2_RADIATION,L_TRACER_MASS,L_CO2_MASS,I_TRACERMASS_START,    &
     &  L_METHANE_LOSS,I_METHANE_TRACERS,L_MCF_LOSS,I_MCF_TRACERNUMBER, &
     &  L_RADON_DECAY,I_RADON_TRACERNUMBER,                             &
     &  L_Q10, L_NEG_TSTAR, L_VEG_FRACS, L_TRIFFID, L_PHENOL,           &
     &  L_TRIF_EQ, L_NRUN_MID_TRIF, L_DISTURB,                          &
     &  CAN_MODEL, PHENOL_PERIOD, TRIFFID_PERIOD,                       &
     &  L_USE_SEASALT_INDIRECT, L_USE_BIOGENIC,                         &
     &  L_USE_SEASALT_DIRECT, L_USE_DUST, L_USE_SULPC_INDIRECT_SW,      &
     &  L_USE_SULPC_INDIRECT_LW, L_USE_SULPC_DIRECT,                    &
     &  L_USE_SULPHATE_AUTOCONV, L_USE_SEASALT_AUTOCONV, L_AUTO_DEBIAS, &
     &  L_USE_SULPHATE_SULPC, L_USE_SEASALT_SULPC,                      &
     &  L_OCFF, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_OCFF_AUTOCONV,        &
     &  L_USE_OCFF_SULPC, L_USE_OCFF_DIRECT, L_USE_OCFF_INDIRECT,       &
     &  L_USE_STOCHEM_CH4, L_USE_STOCHEM_O3,                            &
     &  L_MOD_BARKER_ALBEDO, L_USE_SPEC_SEA, L_MOD_K_FLUX, L_CTILE_FIX, &
     &  L_QUAD_SRC_FIX, L_USE_TPPS_OZONE, I_OZONE_INT,                  &
     &  L_ICOUNT,                                                       &
     &  L_use_orog_corr, L_use_grad_corr,                               &
     &  CALL_CHEM_FREQ, L_USE_SOOT_DIRECT, L_USE_SOOT_INDIRECT,         &
     &  L_USE_SOOT_AUTOCONV, L_USE_SOOT_SULPC, L_USE_BMASS_DIRECT,      &
     &  L_USE_BMASS_INDIRECT, L_USE_BMASS_AUTOCONV, L_USE_BMASS_SULPC,  &
     &  L_USE_ARCLBIOM, L_USE_ARCLBLCK,  L_USE_ARCLSSLT,                &
     &  L_USE_ARCLSULP, L_USE_ARCLDUST,  L_USE_ARCLOCFF, L_USE_ARCLDLTA,&
     &  L_ukca,                                                         &
     &  L_CTILE, L_RIVERS,L_INLAND, RIVER_STEP,                         &
     &  l_sice_meltponds, l_sice_scattering, l_sice_hadgem1a,           &
     &  L_USE_METHOX,Z_TOP,l_mr_physics1,l_mr_physics2,                 &
     &  L_INHOM_CLOUD,                                                  &
     &  L_radiation,L_FORCING,L_TIMESTEP,                               &
     &  L_RADIANCE, L_WENYI, L_bl, L_rain, L_hydrology,                 &
     &  L_TOP,L_PDM,L_USE_AOD,L_USE_CLEARRH,L_UPPER_STRATQ,             &
     &  L_sice_heatflux,L_SOIL_SAT_DOWN,                                &
     &  L_SCVARY,L_VOLCTS,                                              &
     &  l_anthrop_heat_src,                                             &
     &  L_USE_CARIOLLE,L_USE_OZONEINRAD,                                & 
     &  L_RPSEED_READ, L_RPSEED_WRITE,                                  &
     ! Character variables need to be at the end.
     &  H_SECT, A_ASSIM_MODE

      NAMELIST/NLSTCATM/                                                &
     &  Model_domain,L_emcorr,                                          &
     &  L_OASIS, OASIS_COUPLE_FREQ, L_COUPLE_MASTER,                    &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                          &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,              &
     &  problem_number,                                                 &
     &  L_SNOW_ALBEDO,l_ssice_albedo,L_MICROPHY,H_SWBANDS,H_LWBANDS,    &
     &  A_SW_RADSTEP,A_LW_RADSTEP,                                      &
     &  A_SW_RADSTEP_DIAG,A_LW_RADSTEP_DIAG,                            &
     &  A_SW_RADSTEP_PROG,A_LW_RADSTEP_PROG,                            &
     &  L_CLD_AREA,L_ACF_CUSACK,L_ACF_BROOKS,L_PC2,L_PC2_RESET,         &
     &  L_PC2_LBC,L_PC2_diag_sh, L_rad_deg, L_AC, L_VINT_TP,            &
     &  A_ASSIM_START_MIN, A_ASSIM_END_MIN,                             &
     &  NPMSL_HEIGHT,                                                   &
     &  LMASS_CORR , LQT_CORR, LEMQ_PRINT ,A_ENERGYSTEPS,               &
     &  L_GWD,L_USE_USSP,                                               &
     &  L_Murk, L_MURK_ADVECT, L_MURK_SOURCE, L_MURK_BDRY,              &
     &  L_MURK_RAD, L_murk_lbc, L_BL_TRACER_MIX, L_int_uvw_lbc,         &
     &  L_DUST,L_CAM_DUST,L_SULPC_SO2,L_SULPC_DMS,L_SULPC_OZONE,        &
     &  L_SULPC_SO2_O3_NONBUFFERED, L_SO2_SURFEM, L_SO2_HILEM,          &
     &  L_SO2_NATEM,                                                    &
     &  L_SULPC_ONLINE_OXIDANTS,                                        &
     &  L_DMS_EM, L_DMS_EM_INTER,                                       &
     &  L_DMS_Ointer,                                                   &
     &  L_DMS_Liss_Merlivat, L_DMS_Wanninkhof, L_DMS_Nightingale,       &
     &  L_SULPC_NH3, L_NH3_EM, L_SOOT, L_SOOT_SUREM, L_SOOT_HILEM,      &
     &  L_BIOMASS, L_BMASS_SUREM, L_BMASS_HILEM,                        &
     &  L_PSSWND, L_DUST2OCN,                                           &
     &  L_RHCPT, L_CCRAD, L_3D_CCA, L_3D_CCW, L_PHASE_LIM,              &
     &  L_CO2_INTERACTIVE,L_CO2_EMITS,l_cable,                          &
! rml 1/7/13
     &  L_CO2_TRACER,                                                   &
!kdcorbin, 08/10
     &  L_CO2_RADIATION,L_TRACER_MASS,L_CO2_MASS,I_TRACERMASS_START,    &
     &  L_METHANE_LOSS,I_METHANE_TRACERS,L_MCF_LOSS,I_MCF_TRACERNUMBER, &
     &  L_RADON_DECAY,I_RADON_TRACERNUMBER,                             &  
     &  L_Q10, L_NEG_TSTAR, L_VEG_FRACS, L_TRIFFID, L_PHENOL,           &
     &  L_TRIF_EQ, L_NRUN_MID_TRIF, L_DISTURB,                          &
     &  CAN_MODEL, PHENOL_PERIOD, TRIFFID_PERIOD,                       &
     &  L_USE_SEASALT_INDIRECT, L_USE_BIOGENIC,                         &
     &  L_USE_SEASALT_DIRECT, L_USE_DUST, L_USE_SULPC_INDIRECT_SW,      &
     &  L_USE_SULPC_INDIRECT_LW, L_USE_SULPC_DIRECT,                    &
     &  L_USE_SULPHATE_AUTOCONV, L_USE_SEASALT_AUTOCONV, L_AUTO_DEBIAS, &
     &  L_USE_SULPHATE_SULPC, L_USE_SEASALT_SULPC,                      &
     &  L_OCFF, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_OCFF_AUTOCONV,        &
     &  L_USE_OCFF_SULPC, L_USE_OCFF_DIRECT, L_USE_OCFF_INDIRECT,       &
     &  L_USE_STOCHEM_CH4, L_USE_STOCHEM_O3,                            &
     &  L_MOD_BARKER_ALBEDO, L_USE_SPEC_SEA, L_MOD_K_FLUX, L_CTILE_FIX, &
     &  L_QUAD_SRC_FIX, L_USE_TPPS_OZONE, I_OZONE_INT,                  &
     &  L_ICOUNT,                                                       &
     &  L_use_orog_corr, L_use_grad_corr,                               &
     &  CALL_CHEM_FREQ, L_USE_SOOT_DIRECT, L_USE_SOOT_INDIRECT,         &
     &  L_USE_SOOT_AUTOCONV, L_USE_SOOT_SULPC, L_USE_BMASS_DIRECT,      &
     &  L_USE_BMASS_INDIRECT, L_USE_BMASS_AUTOCONV, L_USE_BMASS_SULPC,  &
     &  L_USE_ARCLBIOM, L_USE_ARCLBLCK,  L_USE_ARCLSSLT,                &
     &  L_USE_ARCLSULP, L_USE_ARCLDUST,  L_USE_ARCLOCFF, L_USE_ARCLDLTA,&
     &  L_ukca,                                                         &
     &  L_CTILE, L_RIVERS,L_INLAND, RIVER_STEP,                         &
     &  l_sice_meltponds, l_sice_scattering, l_sice_hadgem1a,           &
     &  L_USE_METHOX,Z_TOP,l_mr_physics1,l_mr_physics2,                 &
     &  L_INHOM_CLOUD,                                                  &
     &  L_radiation,L_FORCING,L_TIMESTEP,                               &
     &  L_RADIANCE, L_WENYI, L_bl, L_rain, L_hydrology,                 &
     &  L_TOP,L_PDM,L_USE_AOD,L_USE_CLEARRH,L_UPPER_STRATQ,             &
     &  L_SCVARY,L_VOLCTS,                                              &
     &  L_sice_heatflux, L_SOIL_SAT_DOWN,                               &
     &  l_anthrop_heat_src,                                             &
     &  L_USE_CARIOLLE,L_USE_OZONEINRAD,                                &
     &  L_RPSEED_READ, L_RPSEED_WRITE,                                  &
     &  H_SECT, A_ASSIM_MODE

      ! Control switches derived and set within the model, hence not
      ! passed in from the UMUI via namelist:

      ! Radiation
      LOGICAL :: L_SW_RADIATE  ! Activate SW radiation this timestep
      LOGICAL :: L_LW_RADIATE  ! Activate LW radiation this timestep
      LOGICAL :: L_SW_RADIATE_DIAG
! Activate fast SW radiation this timestep (3C)
      LOGICAL :: L_LW_RADIATE_DIAG
! Activate fast LW radiation this timestep (3C)
      LOGICAL :: L_SW_RADIATE_PROG
! Activate slow SW radiation this timestep (3C)
      LOGICAL :: L_LW_RADIATE_PROG
! Activate slow LW radiation this timestep (3C)
      LOGICAL :: Lexpand_ozone ! Convert zonal mean ozone to field

      ! convert zonal mean tpps ozone to field
      LOGICAL :: Lexpand_tpps_ozone

      INTEGER :: I_tpps_ozone_opts ! options for tropopause-based ozone

      ! size of super array holding all tracers
      Integer ::   super_array_size
      Integer ::   moisture_array_size

      COMMON/CNTLCATM2/                                                 &
       !  super tracer array size
     &  super_array_size, moisture_array_size,                          &
     &  L_SW_RADIATE,L_LW_RADIATE,                                      &
     &  L_SW_RADIATE_DIAG,L_LW_RADIATE_DIAG,                            &
     &  L_SW_RADIATE_PROG,L_LW_RADIATE_PROG,                            &
     &  Lexpand_ozone,                                                  &
     & Lexpand_tpps_ozone, I_tpps_ozone_opts
! options for tropopause-based ozone



      LOGICAL :: LTLEADS ! Switch for Leads temperature.
!                      ! If FALSE, they are assumed to be TFS
!                      ! Else they are prognostic.
! Default setting: leads temperatures are set to TFS
! HARDWIRE LTLEADS
      PARAMETER(LTLEADS =.FALSE.)
! ----------------------- Comdeck: CNTLOCN  ----------------------------
! Description: COMDECK defining Control variables for the Ocean
!              internal model.
!   This comdeck contains logical variables which are used on the
!   control of certain sections of Ocean model code
!   They replace the previous method of controlling code using *IF DEFs.
!
! Author : R.T.H.Barnes & R.Hill
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER :: O_CLM_START_HR    ! Time ocean climate increments start
      INTEGER :: O_CLM_END_HR      ! Time ocean climate increments end
      INTEGER :: O_INT_CLM_INC     ! # ocean steps  } climate incs.
      INTEGER :: O_INT_ANA_STP     ! # between      } analysis steps

      ! # ocean steps between fwd evolution  of bathys and tesacs
      INTEGER :: O_INT_EVO_BTS

      ! # ocean steps between re-calculation of future bathys and tesacs
      ! valid at this hour

      INTEGER :: O_INT_VRY_BTS

      INTEGER :: O_INT_WTS_ACC    ! # ocean steps betwn accumulating wts

      INTEGER :: O_INT_OBS_FRSH   ! # ocean  } reading new OBS files
      INTEGER :: O_INT_OBS_OUT    ! # steps  } outputting new OBS files
      INTEGER :: O_INT_OBS_STR    ! # between} caching OBS array
      INTEGER :: O_INT_FLD_STR    ! #        } caching model fields

      ! Time at which data assimilation starts (Hours after Basis Time)
      INTEGER :: O_ASSIM_START_HR

      ! Time at which data assimilation ends (Hours after Basis Time)
      INTEGER :: O_ASSIM_END_HR

      INTEGER :: O_ASSIM_ANAL_PER ! Period between analyses (Hours)
      INTEGER :: O_ASSIM_1ST_ANAL    ! First analysis time (Hours)
      LOGICAL :: L_FLUXCORR   ! Heat & water flux correction
      LOGICAL :: L_OGLOBAL    ! Global ocean
      LOGICAL :: L_ICEEVP    ! Use Elastic-Viscous-Plastic Ice dynamics
      LOGICAL :: L_ICESSTILT ! Include sea-surface tilt forcing on ice
      LOGICAL :: L_ICYNPOL   ! North pole fix for use with dynamic ice.
      LOGICAL :: L_ICEFREEDR  ! Free Drift Sea Ice model
      LOGICAL :: L_ICESIMPLE  ! Simple Advection Sea Ice model
      LOGICAL :: L_HADCM4O2I  ! HADCM4 version of ocean-to-ice heat flux
      LOGICAL :: L_IFILTER    ! Filter ice velocities
      LOGICAL :: L_IMCPHEE    ! McPhee ocean-to-ice heat flux
      LOGICAL :: L_IHANEY     ! Haney Forcing Ice
      LOGICAL :: L_ICEITD     ! Use ice thickness distribution
                              ! (i.e. multiple ice categories)
      LOGICAL :: L_ICONCHECK  ! Check ice conservation
                              ! (always FALSE if L_ICEITD=FALSE)
      LOGICAL :: L_ISTRHSM    ! Smooth ice strength
      LOGICAL :: L_ISTRH_PSTAR! If true, use Hibler 79 ice strength
                              ! formula (uses pstar), else use
                              ! Rothrock 75 ice strength formula.
      LOGICAL :: L_ICEFLUXSC  ! Scale A-O coupling of ice related fluxes
                              !    using ice concentrations
      LOGICAL :: L_OADGHR2    ! Ocean assimilation diagnostics
      LOGICAL :: L_OBDY_NORTH   ! Update northern lateral boundary
      LOGICAL :: L_OBDY_SOUTH   ! Update southern lateral boundary
      LOGICAL :: L_OBDY_EAST    ! Update eastern lateral boundary
      LOGICAL :: L_OBDY_WEST    ! Update western lateral boundary
      LOGICAL :: L_OGILL_LBCS   ! Use the Gill boundary scheme
      LOGICAL :: L_OFRS_LBCS    ! Use the FRS boundary scheme
      LOGICAL :: L_OSTVNS_LBCS  ! Use the Stevens boundary scheme
      LOGICAL :: L_OBDY_TRACER  ! Update the tracers
      LOGICAL :: L_OBDY_UV      ! Update the velocities
      LOGICAL :: L_OBDY_STREAM  ! Update the stream functions
      LOGICAL :: L_OBDY_ICE     ! Update ice fields (snow, aice, hice)
!  Start of switches for ocean biogeochemistry and air-sea gas flux
      LOGICAL :: L_OCARBON    ! Carbon cycle model
      ! interactive 3D CO2 field for use with carbon cycle model
      LOGICAL :: L_CO2O_INTERACTIVE
      LOGICAL :: L_OCARB14    ! Calculate atmospheric C12/C14 ratio
      LOGICAL :: L_OCSCEN     ! have a scenario for atmosphere CO2
      LOGICAL :: L_OANCACO2   ! read atmospheric CO2 from ancillary
      LOGICAL :: L_OEXTRAC    ! have two carbon tracers for scenario
      LOGICAL :: L_OCO2OCMIP  ! use carbo chem equm consts from OCMIP
      LOGICAL :: L_OALKSAL    ! base alkalinity on salin (if no bio)
      LOGICAL :: L_OVIRTUAL   ! include virtual surface fluxes
      LOGICAL :: L_OBIOLOGY   ! Effect of phytoplankton on carbon cycle
      LOGICAL :: L_ONUTRELAX  ! relaxation of nutrient to levitus
      LOGICAL :: L_OFEBIO     ! switch for iron-limitation scheme
      LOGICAL :: L_DUST2OCN_O ! get dust deposition from the atmosphere
      LOGICAL :: L_O2BLM4BIO  ! use 2-band light model for bio PrimProd
      LOGICAL :: L_OBIOSTD    ! use standard HadOCC biology
      LOGICAL :: L_OBIODETC   ! separate detrl C,N varbls in std HadOCC
      LOGICAL :: L_OBIODOM    ! use DOM biology
      LOGICAL :: L_OBIODTM    ! switch for Diatom model
      LOGICAL :: L_SWTCHGRZ   ! switch for Fasham switching grazer
      LOGICAL :: L_OFRATIO    ! do ammonium calculations
      LOGICAL :: L_OSRFLX     ! put detritus reaching bottom in surface
      LOGICAL :: L_OBIOTLIM   ! have temperature limitation of phyto.
      LOGICAL :: L_OSHLWCO3   ! even shallow water columns form CaCO3
      LOGICAL :: L_OCTOCHL    ! variable C:Chl for phytoplankton
      LOGICAL :: L_OADVCHL    ! advect chlorophyll as a tracer
      LOGICAL :: L_OCCHLPANC  ! Read carbon:chl ratio from ancillary
      LOGICAL :: L_OOXYGEN    ! Include Oxygen tracer
      LOGICAL :: L_ODMS       ! calculate DMS and flux to the atmosphere
      LOGICAL :: L_OBDIAGS2   ! use alternative bio-diagnostics (DTM)
      LOGICAL :: L_OEXTRACERS ! Include extra tracers
      LOGICAL :: L_OC14TRAC   ! Include Carbon 14 tracer
      LOGICAL :: L_OBOMC14    ! Include bomb Carbon 14 tracer
      LOGICAL :: L_OCFC1112   ! run CFCs as tracers
      LOGICAL :: L_OHELIUM    ! Include Helium-3 and Helium-4 tracers
      LOGICAL :: L_ANCWND_O   ! read 10m windspeed from an ancillary
      LOGICAL :: L_PSSWND_O   ! 10m windspeed passed from atmosphere
      LOGICAL :: L_ICECMSK    ! read an ice mask from an ancillary
      LOGICAL :: L_OCHLANC    ! Read surface chlorophyll from ancillary
      LOGICAL :: L_OLISS      ! Liss & Merlivat wind mixing of tracers
      LOGICAL :: L_OLISS660   ! in Liss/Mer normalise to 660 (old,wrong)
      LOGICAL :: L_OWKHOF     ! Use Wanninkhof 92 piston vel scheme
      LOGICAL :: L_ONGALE     ! Use Nightingale & al piston vel scheme
      LOGICAL :: L_OMNTHWND   ! Use monthly winds in piston vel calc
      LOGICAL :: L_ODLYWND    ! Use daily winds in piston vel calc
!  End of switches for ocean biogeochemistry and air-sea gas flux
      LOGICAL :: L_OCNASSM    ! Activate ocean assimilation
      LOGICAL :: L_OCYCLIC    ! Cyclic boundary conditions
      LOGICAL :: L_OFILTER    ! Fourier filtering for high latitudes
      LOGICAL :: L_OFILTHARD ! Extra stringency on F. filtering
      LOGICAL :: L_OFILTTROP  ! F. filtering on FS barotropic velocities
      LOGICAL :: L_OFILTLBAL  ! Control F. filter load balancing
      LOGICAL :: L_OFREESFC   ! Use free surface conditions
      LOGICAL :: L_OFSMARGINAL   ! Control IFS marginal seas height
      LOGICAL :: L_FLUXD
      LOGICAL :: L_OHANEY     ! Haney Forcing heat/fresh water fluxes
      LOGICAL :: L_OHMEAD     ! Mead tracer transport diagnostics
      LOGICAL :: L_OICECOUP   ! Coupled model with Sea Ice
      LOGICAL :: L_OFLXNRM    ! Flux inputs normalised over sea-ice
      LOGICAL :: L_OIMPDIF    ! CN vertical diffusion scheme
      LOGICAL :: L_OISLANDS   ! Include Island Routines
      LOGICAL :: L_OISOPYC    ! Isopycnal diffusion scheme
      LOGICAL :: L_OLATVISC   ! Latitude dependent viscosity
      LOGICAL :: L_OANIVISC   ! Anisotropic viscosity (as in GloSea)
      LOGICAL :: L_OMIXLAY    ! Wind mixing of tracers-mixed layer scheme
      LOGICAL :: L_ONOCLIN    ! Barotropic solution
      LOGICAL :: L_ONOPOLO    ! No sea ice at North Pole
      LOGICAL :: L_OPENBC     ! Read in lateral boundary fields
      LOGICAL :: L_ORICHARD   ! Evaluate & use Richardson No.
      LOGICAL :: L_OROTATE    ! Coriolis force calculation
      LOGICAL :: L_OSOLAR     ! Calc solar penetration for given water ty
      LOGICAL :: L_OSOLARAL   ! Calc sol. pen. - simplified layer structu
      LOGICAL :: L_OSYMM      ! Symmetric boundary conditions
      LOGICAL :: L_OVARYT     ! Varying time step with depth
      LOGICAL :: L_ORIVERS    ! River run-off routines
      LOGICAL :: L_SEAICE     ! Include Sea Ice model
      LOGICAL :: L_TRANGRID   ! Spatial interp. in coupled model
      LOGICAL :: L_OCONJ     ! Whether to use conjugate gradient solver
      LOGICAL :: L_UPWIND     ! Upwind differencing for tracer advection
      LOGICAL :: L_OPRINT     ! Whether to print incidental ocean info
      LOGICAL :: L_ODELPLUS   !
      LOGICAL :: L_OTROPIC    !
      LOGICAL :: L_OISOMOM
      LOGICAL :: L_OISOGMSKEW
      LOGICAL :: L_OISOGM
      LOGICAL :: L_OBIHARMGM
      LOGICAL :: L_OBIGMCUBE  ! Cubic cos(lat) term for biharm GM
      LOGICAL :: L_OVISHADGEM1
      ! Mediterranean outflow - 288*144 and 96*73 grids only - uses
      ! hardwired gridpoint nos
      LOGICAL :: L_OMEDOUT
      LOGICAL :: L_OCONVROUS  ! Roussenov convective adjustment
      LOGICAL :: L_OEXTRAP ! Extrapolation of vertical density gradients
      LOGICAL :: L_OISOPYCGM  ! Gent and McWilliams eddy parametrisation
      LOGICAL :: L_OISOTAPER  ! Tapering of isopycnal diffusion
      LOGICAL :: L_OVISBECK    ! Visbeck scheme
      LOGICAL :: L_OVISPLUS    ! Enhanced Visbeck for high lat damping
      LOGICAL :: L_OQLARGE     ! Quadratic Large scheme
      LOGICAL :: L_OFULARGE   ! FULL LARGE SCHEME
      LOGICAL :: L_OPANDP     ! RI-DEPENDENT VERT MIX SCHEMES
      LOGICAL :: L_OSTATEC    ! DENSITY CHOICE FOR RI-CALC
      LOGICAL :: L_OUSTARWME  ! WME OR WSTRESS TO FIND USTAR

      LOGICAL :: L_OZVRT      ! barotropic vorticity diagnostic switch
                           ! set by OCN_FOR_STEP (not in namelist)
      LOGICAL :: L_SLOPEMAX   ! Selects SLOPE_MAX isopycnal diffusion
      LOGICAL :: L_COXCNVC    ! Selects original Cox convection scheme
      LOGICAL :: L_OMEDADV
      LOGICAL :: L_OHUDOUT
      LOGICAL :: L_OSTRAIT    ! T => Strait exchange flows parametrised
      LOGICAL :: L_OSTR_CLM   ! T => some strait pts set by climate
                              !      values
      LOGICAL :: L_REFSAL
      LOGICAL :: L_SALFLUXFIX
      LOGICAL :: L_INLANSEA
      LOGICAL :: L_OBOTFRI
      LOGICAL :: L_OEOS25
      LOGICAL :: L_OBIMOM  ! biharmonic momentum diffusion
      LOGICAL :: L_OBISURF  ! biharmonic tracer diffusion in top layers
      LOGICAL :: L_OBMCUBE ! Cubic cos(lat) term for biharm mom diff
      LOGICAL :: L_OBULKRI
      LOGICAL :: L_OWINDMIX
      LOGICAL :: L_OBULKMAXMLD
      LOGICAL :: L_OBDBBL
      LOGICAL :: L_OBIAS
      LOGICAL :: L_ORLP       ! Select rigid lid pressure calculation
      ! Additions to CCONTROL for ocean assimilation

      LOGICAL :: LAS_CLM_INC   ! make increments to relax to climate
      LOGICAL :: LAS_ADD_INC   ! add or subtract analysis increments
      LOGICAL :: LAS_ANA_STP   ! calculate analysis increments
      LOGICAL :: LAS_EVO_BTS   ! evolve bathy and tesac obs 1 step
      LOGICAL :: LAS_VRY_BTS   ! estimate bathys and tesacs at this hour
      LOGICAL :: LAS_WTS_ACC   ! evolve accumulated weights
      LOGICAL :: LAS_OBS_FRSH  ! to refresh main OBS data set
      LOGICAL :: LAS_OBS_OUT   ! output ACOBS file for incremented obs
      LOGICAL :: LAS_FLD_STR   ! output model fields to cache store
      LOGICAL :: LAS_OBS_STR   ! output obs to cache store
      LOGICAL :: L_OMNRLP      ! =T if mean rlp to be read from dump
      LOGICAL :: L_OHADGEM1      ! controls HADGEM1 specific code


      NAMELIST / NLSTCOCN /                                             &
     &  O_CLM_START_HR, O_CLM_END_HR, O_INT_CLM_INC, O_INT_ANA_STP,     &
     &  O_INT_EVO_BTS, O_INT_VRY_BTS, O_INT_WTS_ACC, O_INT_OBS_FRSH,    &
     &  O_INT_OBS_OUT, O_INT_OBS_STR, O_INT_FLD_STR,                    &
     &  O_ASSIM_START_HR, O_ASSIM_END_HR, O_ASSIM_ANAL_PER,             &
     &  O_ASSIM_1ST_ANAL,L_FLUXCORR,L_OGLOBAL,L_ISTRHSM,L_ISTRH_PSTAR,  &
     &  L_IFILTER,L_ICEFLUXSC,                                          &
     &  L_ICEEVP,L_ICESSTILT,L_ICYNPOL,L_ICEFREEDR, L_ICESIMPLE,        &
     &  L_IHANEY, L_ICEITD, L_ICONCHECK,L_HADCM4O2I,L_IMCPHEE,          &
     &  L_OADGHR2,L_OBDY_NORTH,L_OBDY_SOUTH,L_OBDY_EAST,L_OBDY_WEST,    &
     &  L_OGILL_LBCS,L_OFRS_LBCS,L_OSTVNS_LBCS,                         &
     &  L_OBDY_TRACER,L_OBDY_UV,L_OBDY_STREAM,L_OBDY_ICE,               &
!  Start of switches for ocean biogeochemistry and air-sea gas flux
     &  L_OCARBON, L_CO2O_INTERACTIVE, L_OCARB14, L_OCSCEN, L_OANCACO2, &
     &  L_OEXTRAC, L_OCO2OCMIP, L_OALKSAL, L_OVIRTUAL,                  &
     &  L_OBIOLOGY, L_ONUTRELAX, L_OFEBIO, L_DUST2OCN_O, L_O2BLM4BIO,   &
     &  L_OBIOSTD, L_OBIODETC, L_OBIODOM, L_OBIODTM, L_SWTCHGRZ,        &
     &  L_OFRATIO, L_OSRFLX, L_OBIOTLIM, L_OSHLWCO3,                    &
     &  L_OCTOCHL, L_OADVCHL, L_OCCHLPANC,                              &
     &  L_OOXYGEN, L_ODMS, L_OBDIAGS2,                                  &
     &  L_OEXTRACERS, L_OC14TRAC, L_OBOMC14, L_OCFC1112, L_OHELIUM,     &
     &  L_ANCWND_O, L_PSSWND_O, L_ICECMSK, L_OCHLANC,                   &
     &  L_OLISS, L_OLISS660, L_OWKHOF, L_ONGALE, L_OMNTHWND, L_ODLYWND, &
!  End of switches for ocean biogeochemistry and air-sea gas flux
     &  L_OCNASSM, L_OCYCLIC, L_OFILTER, L_OFREESFC, L_OFSMARGINAL,     &
     &  L_FLUXD, L_OFILTLBAL, L_OFILTHARD, L_OFILTTROP,                 &
     &  L_OHANEY, L_OHMEAD, L_OICECOUP, L_OFLXNRM,                      &
     &  L_OIMPDIF, L_OISLANDS, L_OISOPYC, L_OLATVISC, L_OANIVISC,       &
     &  L_OMIXLAY, L_ONOCLIN, L_ONOPOLO, L_OPENBC,                      &
     &  L_ORICHARD, L_OROTATE, L_OSOLAR, L_OSOLARAL,                    &
     &  L_OSYMM, L_OVARYT, L_ORIVERS, L_SEAICE, L_OCONJ,                &
     &  L_TRANGRID, L_UPWIND, L_OPRINT,                                 &
     &  L_ODELPLUS, L_OTROPIC,                                          &
     &  L_OBIGMCUBE,                                                    &
     &  L_OISOMOM,L_OISOGMSKEW,L_OISOGM,L_OBIHARMGM,L_OVISHADGEM1,      &
     &  L_OMEDOUT,                                                      &
     &  L_OCONVROUS,                                                    &
     &  L_OEXTRAP,L_OISOPYCGM,L_OISOTAPER,                              &
     &  L_OVISBECK,                                                     &
     &  L_OVISPLUS,                                                     &
     & L_OBIMOM, L_OBISURF,                                             &
     &  L_OBMCUBE,                                                      &
     &  L_OQLARGE,                                                      &
     &  L_OMEDADV,L_OHUDOUT, L_OSTRAIT, L_OSTR_CLM,                     &
     &  L_REFSAL,L_SALFLUXFIX,L_INLANSEA,L_OBOTFRI,                     &
     &  L_OEOS25,                                                       &
     &  L_OFULARGE,L_OPANDP,L_OSTATEC,L_OUSTARWME,                      &
     &  L_SLOPEMAX,L_COXCNVC,                                           &
     &  L_OBDBBL,                                                       &
     &  L_OBIAS,                                                        &
     &  L_ORLP,                                                         &
      ! additions for control of ocean assimilation
     &  LAS_ADD_INC,LAS_CLM_INC,LAS_ANA_STP,                            &
     &  LAS_EVO_BTS,LAS_VRY_BTS,LAS_WTS_ACC,                            &
     &  LAS_OBS_FRSH,LAS_OBS_OUT,LAS_FLD_STR,LAS_OBS_STR,               &
     &  L_OMNRLP, L_OHADGEM1

      COMMON / CNTLCOCN /                                               &
     &  O_CLM_START_HR, O_CLM_END_HR, O_INT_CLM_INC, O_INT_ANA_STP,     &
     &  O_INT_EVO_BTS, O_INT_VRY_BTS, O_INT_WTS_ACC, O_INT_OBS_FRSH,    &
     &  O_INT_OBS_OUT, O_INT_OBS_STR, O_INT_FLD_STR,                    &
     &  O_ASSIM_START_HR, O_ASSIM_END_HR, O_ASSIM_ANAL_PER,             &
     &  O_ASSIM_1ST_ANAL,L_FLUXCORR,L_OGLOBAL,L_ISTRHSM,L_ISTRH_PSTAR,  &
     &  L_IFILTER,L_ICEFLUXSC,                                          &

     &  L_ICEEVP,L_ICESSTILT,L_ICYNPOL,L_ICEFREEDR, L_ICESIMPLE,        &
     &  L_IHANEY, L_ICEITD, L_ICONCHECK, L_HADCM4O2I,L_IMCPHEE,         &
     &  L_OADGHR2,L_OBDY_NORTH,L_OBDY_SOUTH,L_OBDY_EAST,L_OBDY_WEST,    &
     &  L_OGILL_LBCS,L_OFRS_LBCS,L_OSTVNS_LBCS,                         &
     &  L_OBDY_TRACER,L_OBDY_UV,L_OBDY_STREAM,L_OBDY_ICE,               &
!  Start of switches for ocean biogeochemistry and air-sea gas flux
     &  L_OCARBON, L_CO2O_INTERACTIVE, L_OCARB14, L_OCSCEN, L_OANCACO2, &
     &  L_OEXTRAC, L_OCO2OCMIP, L_OALKSAL, L_OVIRTUAL,                  &
     &  L_OBIOLOGY, L_ONUTRELAX, L_OFEBIO, L_DUST2OCN_O, L_O2BLM4BIO,   &
     &  L_OBIOSTD, L_OBIODETC, L_OBIODOM, L_OBIODTM, L_SWTCHGRZ,        &
     &  L_OFRATIO, L_OSRFLX, L_OBIOTLIM, L_OSHLWCO3,                    &
     &  L_OCTOCHL, L_OADVCHL, L_OCCHLPANC,                              &
     &  L_OOXYGEN, L_ODMS, L_OBDIAGS2,                                  &
     &  L_OEXTRACERS, L_OC14TRAC, L_OBOMC14, L_OCFC1112, L_OHELIUM,     &
     &  L_ANCWND_O, L_PSSWND_O, L_ICECMSK, L_OCHLANC,                   &
     &  L_OLISS, L_OLISS660, L_OWKHOF, L_ONGALE, L_OMNTHWND, L_ODLYWND, &
!  End of switches for ocean biogeochemistry and air-sea gas flux
     &  L_OCNASSM, L_OCYCLIC, L_OFILTER, L_OFREESFC, L_OFSMARGINAL,     &
     &  L_FLUXD, L_OFILTLBAL, L_OFILTHARD, L_OFILTTROP,                 &
     &  L_OHANEY, L_OHMEAD, L_OICECOUP, L_OFLXNRM,                      &
     &  L_OIMPDIF, L_OISLANDS, L_OISOPYC, L_OLATVISC, L_OANIVISC,       &
     &  L_OMIXLAY, L_ONOCLIN, L_ONOPOLO, L_OPENBC,                      &
     &  L_ORICHARD, L_OROTATE, L_OSOLAR, L_OSOLARAL,                    &
     &  L_OSYMM, L_OVARYT, L_ORIVERS, L_SEAICE, L_OCONJ,                &
     &  L_TRANGRID, L_UPWIND, L_OPRINT,                                 &
     &  L_ODELPLUS, L_OTROPIC, L_OZVRT,                                 &
     &  L_OMEDOUT,L_OISOMOM,L_OISOGMSKEW,L_OISOGM,                      &
     &  L_OBIGMCUBE,                                                    &
     &  L_OBIHARMGM,L_OVISHADGEM1,                                      &
     &  L_OCONVROUS,                                                    &
     &  L_OEXTRAP,L_OISOPYCGM,L_OISOTAPER,                              &
     &  L_OVISBECK,                                                     &
     &  L_OVISPLUS,                                                     &
     & L_OBIMOM, L_OBISURF,                                             &
     &  L_OBMCUBE,                                                      &
     &  L_OQLARGE,                                                      &
     &  L_OMEDADV,L_OHUDOUT, L_OSTRAIT, L_OSTR_CLM,                     &
     &  L_REFSAL,L_SALFLUXFIX,L_INLANSEA,L_OBOTFRI,                     &
     &  L_OEOS25,                                                       &
     &  L_OFULARGE,L_OPANDP,L_OSTATEC,L_OUSTARWME,                      &
     &  L_OBULKRI,L_OWINDMIX,L_OBULKMAXMLD,                             &
     &  L_SLOPEMAX,L_COXCNVC,                                           &
     &  L_OBDBBL,                                                       &
     &  L_OBIAS,                                                        &
     &  L_ORLP,                                                         &
      ! additions for control of ocean assimilation
     &  LAS_ADD_INC,LAS_CLM_INC,LAS_ANA_STP,                            &
     &  LAS_EVO_BTS,LAS_VRY_BTS,LAS_WTS_ACC,                            &
     &  LAS_OBS_FRSH,LAS_OBS_OUT,LAS_FLD_STR,LAS_OBS_STR,               &
     &  L_OMNRLP, L_OHADGEM1

! ----------------------- Header file CRUNTIMC  -----------------------
! Description: Run-time constants for the Atmosphere model (read only).
!              Contains variables that define parametrization values
!              chosen for atmosphere physics and dynamics schemes.
!              [Note that CNTLATM holds accompanying control switches
!              needed for addressing.]
!
! Current Code Owner: R. Rawlins
!
!------------------   Physics:   --------------------------------------
! Generalised physics switches:


! Boundary layer:

      LOGICAL :: l_use_bl_diag_term
      LOGICAL :: L_SBLeq
      LOGICAL :: L_SBLco
      LOGICAL :: L_LAMBDAM2    ! LambdaM=2*LambdaH (operational setting)
      LOGICAL :: L_FULL_LAMBDAS! Lambdas NOT reduced above NTML_LOCAL+1
      LOGICAL :: L_us_blsol    ! Switch for stable and non-oscillatory
                               !   BL vertical diffusion scheme
      REAL :: alpha_Cd (max_number_alpha_cds)
      REAL :: Muw_SBL
      REAL :: Mwt_SBL
      REAL :: Charnock
      REAL :: Z0HSEA           ! Roughness lengths for heat and moisture
      REAL :: Z0MIZ            !   over sea, marginal ice zone and
      REAL :: Z0SICE           !   sea ice (m)
      REAL :: ALBSNC_NVG(NNVG) ! Snow-covered albedo.
      REAL :: ALBSNF_NVG(NNVG) ! Snow-free albedo.
      REAL :: CATCH_NVG(NNVG)  ! Canopy capacity (kg/m2).
      REAL :: GS_NVG(NNVG)     ! Surface conductance (m/s).
      REAL :: INFIL_NVG(NNVG)  ! Infiltration enhancement factor.
      REAL :: ROOTD_NVG(NNVG)  ! Rootdepth (m).
      REAL :: Z0_NVG(NNVG)     ! Roughness length (m).
      REAL :: CH_NVG(NNVG)     ! "Canopy" heat capacity (J/K/m2)
      REAL :: VF_NVG(NNVG)     ! Fractional "canopy" coverage
      REAL :: Puns, Pstb ! parameters for uncond stable numerical solver
                         ! Puns : used in an unstable BL column
                         ! Pstb : used in an stable BL column
      REAL :: OROG_DRAG_PARAM  ! Drag coefficient for orographic form drag
      LOGICAL :: L_EMIS_LAND_GEN ! Aggregate land surface emissivity
      REAL :: EMIS_LAND_GEN      ! Aggregate land surface emissivity
      REAL :: WeightLouisToLong_in
!                        ! Weighting of Louis tails towards long
!                        ! tails: this should be removed in favour
!                        ! of placing RUN_BL in the module bl_option_mod 
!                        ! at a later release.

      INTEGER, DIMENSION(20) :: BL_OPTIONS
      INTEGER :: DECFIX
      INTEGER :: STOPWE_SBL
      INTEGER :: TRWEIGHTS1
      INTEGER :: FLUX_GRAD
      INTEGER :: NG_STRESS
      INTEGER :: SBL_OP
      INTEGER :: ISHEAR_BL
      INTEGER :: FORMDRAG, COR_UST, COR_MO_ITER
      INTEGER :: NON_LOCAL_BL
      INTEGER :: LOCAL_FA
      INTEGER :: PRANDTL
      INTEGER :: I_SCRN_T_DIAG
      INTEGER :: Buddy_sea
      INTEGER :: FD_stab_dep
      INTEGER :: Keep_Ri_FA
      INTEGER :: NL_BL_LEVELS

      INTEGER :: CAN_RAD_MOD              ! For radiation model (Pete Falloon)
      INTEGER :: ILAYERS                  ! For radiation model (Pete Falloon)

      COMMON  /RUN_BLVEG/ALBSNC_NVG,ALBSNF_NVG,CATCH_NVG,GS_NVG,        &
     &  INFIL_NVG,ROOTD_NVG,Z0_NVG,CH_NVG,VF_NVG,                       &
     &  CAN_RAD_MOD,ILAYERS
       NAMELIST/RUN_BLVEG/ALBSNC_NVG,ALBSNF_NVG,CATCH_NVG,GS_NVG,       &
     &  INFIL_NVG,ROOTD_NVG,Z0_NVG,CH_NVG,VF_NVG,                       &
     &  CAN_RAD_MOD,ILAYERS

      COMMON  /RUN_BLICE/Z0HSEA,Z0MIZ,Z0SICE
      NAMELIST /RUN_BLICE/Z0HSEA,Z0MIZ,Z0SICE
      REAL :: SeaSalinityFactor
!       Scaling of qsat allowing for salinity of sea water
      INTEGER :: ISeaZ0T
!       Option for specifying the thermal roughness length
!       at sea points
      INTEGER :: ISeaDynDiag
!       Option for dynamic diagnosis of boundary layer types
       ! Maximum and minimum values for the STPH_RP scheme
       ! Boundary Layer
       REAL         :: par_mezcla_max,par_mezcla,par_mezcla_min ! Max, 
                      ! mean and min value for the neutral mixing length
       REAL         :: G0_max,G0_RP,G0_min ! Max,mean and min values
                      ! for the flux profile parameter
       REAL         :: Charnock_max,Charnock_min ! Max and min values for
                      ! the charnock parameter

      COMMON  /RUN_BL/ FORMDRAG, OROG_DRAG_PARAM, l_use_bl_diag_term         &
     &  ,L_SBLeq, L_SBLco, alpha_Cd, Muw_SBL, Mwt_SBL, BL_OPTIONS            &
     &  ,par_mezcla_max, par_mezcla, par_mezcla_min                          &
     &  ,G0_max, G0_RP, G0_min, Charnock_max, Charnock, Charnock_min         &
     &  ,L_LAMBDAM2, L_FULL_LAMBDAS, SeaSalinityFactor, ISeaZ0T              &
     &  ,ISeaDynDiag ,L_us_blsol, Puns, Pstb, L_EMIS_LAND_GEN, EMIS_LAND_GEN &
     &  ,WeightLouisToLong_in

      NAMELIST/RUN_BL/ FORMDRAG, OROG_DRAG_PARAM                             &
     &  ,l_use_bl_diag_term, SBL_OP, COR_UST,  COR_MO_ITER                   &
     &  ,L_SBLeq, L_SBLco, alpha_Cd, Muw_SBL, Mwt_SBL, ISHEAR_BL, NG_STRESS  &
     &  ,DECFIX, STOPWE_SBL, TRWEIGHTS1, FLUX_GRAD                           &
     &  ,NON_LOCAL_BL, NL_BL_LEVELS                                          &
     &  ,par_mezcla_max, par_mezcla, par_mezcla_min                          &
     &  ,G0_max, G0_RP, G0_min, Charnock_max, Charnock, Charnock_min         &
     &  ,L_LAMBDAM2, L_FULL_LAMBDAS, LOCAL_FA, Prandtl, SeaSalinityFactor    &
     &  ,ISeaZ0T, ISeaDynDiag, Buddy_sea, L_us_blsol, Puns, Pstb             &
     &  ,FD_stab_dep, Keep_Ri_FA                                             &
     &  ,L_EMIS_LAND_GEN, EMIS_LAND_GEN, I_SCRN_T_DIAG                       &
     &  ,WeightLouisToLong_in

! Surface parameters for plant functional types

      REAL    :: ALBSNC_MAX(NPFT)
      REAL    :: ALBSNC_MIN(NPFT)
      REAL    :: ALBSNF_MAX(NPFT)
      REAL    :: DZ0V_DH(NPFT)
      REAL    :: CATCH0(NPFT)
      REAL    :: DCATCH_DLAI(NPFT)
      REAL    :: INFIL_F(NPFT)
      REAL    :: KEXT(NPFT)
      REAL    :: ROOTD_FT(NPFT)

      COMMON  /RUN_PFT/ALBSNC_MAX,ALBSNC_MIN,ALBSNF_MAX,DZ0V_DH,        &
     &  CATCH0,DCATCH_DLAI,INFIL_F,KEXT,ROOTD_FT
      NAMELIST/RUN_PFT/ALBSNC_MAX,ALBSNC_MIN,ALBSNF_MAX,DZ0V_DH,        &
     &  CATCH0,DCATCH_DLAI,INFIL_F,KEXT,ROOTD_FT


      ! Large scale precipitation:

      ! Threshold cloud liquid water content over land/sea for
      ! conversion to ppn (kg water per m**3)

      ! Defunct: not supported
      REAL,PARAMETER:: CW_LAND = 8.0E-4
      REAL,PARAMETER:: CW_SEA  = 2.0E-4


      LOGICAL :: L_cry_agg_dep ! Limit deposition to moisture available
      LOGICAL :: L_it_melting  ! Use Iterative melting
      LOGICAL :: L_seq_mcr     ! Use sequential microphysics updating
      LOGICAL :: L_autoc_3b    ! Use 3B autoconversion method
      LOGICAL :: L_autolim_3b  ! Use 3B autoconversion limit
      LOGICAL :: L_autoconv_murk  ! Use murk aerosol to calculate
                                  ! droplet number
      LOGICAL :: L_droplet_settle ! Allow cloud droplets to settle
      LOGICAL :: L_psd         ! Use generic ice particle size distn.

      REAL :: ec_auto          ! Collision / collection efficiency
      REAL :: N_drop_land      ! Droplet concentration over land
      REAL :: N_drop_sea       ! Droplet concentration over sea
      REAL :: Ntot_land        ! Droplet concentration over land
      REAL :: Ntot_sea         ! Droplet concentration over sea
      REAL :: N_drop_land_cr   ! N_drop_land ^ (-1/3)
      REAL :: N_drop_sea_cr    ! N_drop_sea ^ (-1/3)
      REAL :: X1R              ! Intercept parameter for raindrop
                               ! size distribution
      REAL :: X2R              ! Scaling parameter for rain DSD.
      REAL :: X4R              ! Shape factor for rain DSD.
      REAL :: X1I              ! Intercept parameter for aggregate
                               ! size distribution
      REAL :: X1IC             ! Intercept parameter for crystal
                               ! size distribution
      REAL :: AI,  BI          ! Aggregate mass-size params m(D)=AI D^BI
      REAL :: AIC, BIC         ! Crystal mass-size params m(D)=AIC D^BIC
      REAL :: LSP_EI, LSP_FI   ! Aggregate Best-Reynolds relationship
                               ! Re = LSP_EI BE^LSP_FI
      REAL :: LSP_EIC, LSP_FIC ! Crystal Best-Reynolds relationship
                               ! Re = LSP_EIC BE^LSP_FIC

      ! Critical RH for layer cloud formation
      REAL :: RHCRIT(max_model_levels)
       ! Maximum and minimum values for the STPH_RP scheme
       ! Large Scale Precipitation
      REAL         :: RHCRIT_max,RHCRIT_min ! Max and min values
                       ! for the critical relative humidity
      REAL         :: CI_max,CI_min ! Max and min values for the
                       ! CI parameter, which controls the ice fall speed
      REAL         :: M_CI_max,M_CI,M_CI_min ! Max, mean, min values for
                       ! the multiplication factor for CI and CIC in 3C
                       ! Microphysics
 
      COMMON  /RUN_Precip/RHCRIT                                        &
     & ,l_cry_agg_dep,L_it_melting,L_seq_mcr,L_autoc_3b,L_autolim_3b,   &
     &  l_psd,                                                          &
     &  l_autoconv_murk,ec_auto,N_drop_land,N_drop_sea,Ntot_land,       &
     &  Ntot_sea,N_drop_land_cr,N_drop_sea_cr,X1R,X2R,X4R,X1I,X1IC      &
     &  ,AI,BI,AIC,BIC,LSP_EI,LSP_FI,LSP_EIC,LSP_FIC                    &
     &  ,RHCRIT_max,RHCRIT_min                                          &
     &  ,CI_max,CI_min,M_CI_max,M_CI,M_CI_min,L_droplet_settle
      NAMELIST/RUN_Precip/RHCRIT                                        &
     & ,l_cry_agg_dep,L_it_melting,L_seq_mcr,L_autoc_3b,L_autolim_3b,   &
     &  l_psd,                                                          &
     &  l_autoconv_murk,ec_auto,N_drop_land,N_drop_sea,Ntot_land,       &
     &  Ntot_sea,N_drop_land_cr,N_drop_sea_cr,X1R,X2R,X4R,X1I,X1IC      &
     &  ,AI,BI,AIC,BIC,LSP_EI,LSP_FI,LSP_EIC,LSP_FIC                    &
     &  ,RHCRIT_max,RHCRIT_min                                          &
     &  ,CI_max,CI_min,M_CI_max,M_CI,M_CI_min,L_droplet_settle

      ! Large scale cloud:

      LOGICAL :: L_eacf                ! Use empirically adjusted
                                       ! cloud fraction

      INTEGER :: cloud_fraction_method ! Selects total cloud fraction
                                       ! calculation method
      INTEGER :: ice_fraction_method   ! Selects ice cloud fraction
                                       ! calculation method

      REAL    :: overlap_ice_liquid    ! Generic overlap parameter
                                       ! between ice and liquid phases
      REAL    :: ctt_weight            ! Cloud top temperature weight
      REAL    :: t_weight              ! Local temperature weight
      REAL    :: qsat_fixed            ! Prescribed qsat value
      REAL    :: sub_cld               ! Scaling factor
      REAL    :: dbsdtbs_turb_0        ! PC2 erosion rate / s-1

      COMMON  /RUN_Cloud/L_eacf,cloud_fraction_method,                  &
     &  overlap_ice_liquid,ice_fraction_method,ctt_weight,t_weight,     &
     &  qsat_fixed,sub_cld,dbsdtbs_turb_0

      NAMELIST/RUN_Cloud/L_eacf,cloud_fraction_method,                  &
     &  overlap_ice_liquid,ice_fraction_method,ctt_weight,t_weight,     &
     &  qsat_fixed,sub_cld,dbsdtbs_turb_0


      ! Radiation:

      ! True if climatological aerosol is included.
      LOGICAL :: L_climat_aerosol


      ! True to use real boundary layer heights to specify the boundary
      ! layer aerosol.
      LOGICAL :: L_clim_aero_hgt
      LOGICAL :: L_HadGEM1_Clim_Aero
!             Flag to use HadGEM1 setting for climatological aerosols

      LOGICAL :: L_SEC_VAR ! true if using time varying astronomy
      LOGICAL :: L_EqT     ! True if including the equation of time

      INTEGER :: A_SW_SEGMENTS ! No of batches used in shortwave code
      INTEGER :: A_SW_SEG_SIZE ! Size of sw batches. 
      INTEGER :: A_LW_SEGMENTS ! No of batches used in longwave code
      INTEGER :: A_LW_SEG_SIZE ! Size of lw batches.
      INTEGER :: AERO_BL_LEVELS !Common number of layers taken to be 
!                                occupied by the boundary-layer
!                                aerosol if the boundary layer
!                                depth is not used to determine the 
!                                number separately at each grid-point
!                                In previous versions of the code,
!                                this was taken to be BL_LEVELS 

      Logical :: l_ovrlap    ! Requires l_ccrad=.TRUE.
                             ! Allows Convective and LS Cloud to overlap 
                             ! for radiative impacts.
                             ! (Experimental, defaulted to FALSE,
                             ! requires a hand-edit to change)

      Logical :: l_ccw_scav  ! Requires l_ccrad=.TRUE. .AND. l_ovrlap=.TRUE.
                             ! Allows Convective Cloud Water (CCW) to
                             ! compensate for LS Cloud water in overlapping
                             ! LS/CCA fractions.
                             ! (Experimental, defaulted to FALSE, requires a
                             ! hand-edit to change)

      Logical :: l_cldtop_t_fix
                              ! .TRUE. Corrects Cloud Top Temperature in
                              !        radiation

      ! Valid options are specified in the include file o3intp.h.

      ! radiation block is already a separate module for some variables:
!
! Description:
!   Common Block for radiation
!   Contains:
!   Parameters used to increase the true surface albedo linearly from
!   a minimum value, alphaM at melting point, to a "cold-ice" value,
!   alphaC, over a temperature range, dTice
!   Mass Mixing Ratios of minor Gases N2O,CH4,CFC11,CFC12,O2
!   CFC113,HCFC22,HFC125, and HFC134A
!
! Current Code Owner: S Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.4      09/09/94 Original code. C Hewitt
! 3.4      02/10/94 Taken over, and MMRs added  S Woodward
! 4.4      18/09/97 Mass Mixing Ratio of Oxygen added
! 4.5      21/08/98 ALPHAB introduced and ALPHAM,ALPHAC reinterpreted
!          if l_ssice_albedo i.e. if snow on sea-ice affects albedo.
!          Jonathan Gregory
! 4.5      18/05/98 Mass Mixing Ratios of further halocarbons added.
!                                    J. M. Edwards
! 5.5      11/02/03 Added IS_NCOL    K.D.Williams
! 6.2      06/01/06 sea ice albedo parameters added  J.Ridley
!
! Declarations:
! Imported global variables (*CALLed COMDECKs etc...)

! Global parameters:
! Albedo of sea-ice at melting point (TM) if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at melting point (TM) if l_ssice_albedo
      REAL ALPHAM  ! "M" for "melting"
! Albedo of sea-ice at and below TM-DTICE if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at and below TM-DTICE if l_ssice_albedo
      REAL ALPHAC  ! "C" for "cold"
! Albedo of snow-free sea-ice if l_ssice_albedo
      REAL ALPHAB  ! "B" for "bare"
! Temperature range in which albedo of sea-ice, if .not.l_ssice_albedo,
! or of snow on sea-ice, if l_ssice_albedo, varies between its limits
      REAL DTICE
! Temperature range below TM over which meltponds form if
! l_sice_meltponds and l_ssice_albedo
      REAL DT_BARE
! Increment to albedo for each degree temperature rises above
! TM-DT_BARE. Only used if l_sice_meltponds and l_ssice_albedo
      REAL DALB_BARE_WET
! Fraction of SW radiation that penetrates seaice and scatters
! back causing an addition to the albedo. Only active if l_ssice_albedo
! and l_sice_scattering
      REAL PEN_RAD_FRAC
! attenutation parameter for SW in seaice which controls the
! additional albedo due to internal scattering
      REAL SW_BETA
      REAL N2OMMR  ! N2O mmr (Mass Mixing Ratio)
      REAL CH4MMR  ! CH4 mmr
      REAL C11MMR  ! CFC11 mmr
      REAL C12MMR  ! CFC12 mmr
      REAL O2MMR   ! O2 mmr
      REAL C113MMR      ! CFC113 mmr
      REAL HCFC22MMR    ! HCFC22 mmr
      REAL HFC125MMR    ! HFC125 mmr
      REAL HFC134AMMR   ! HFC134A mmr
! Number of columns used for internal sampling by the ISCCP simulator
      INTEGER IS_NCOL

! Global scalars:

! Global dynamic arrays:

! COMMON blocks:
      COMMON /RAD_COM/ ALPHAM, ALPHAC, ALPHAB, DTICE, N2OMMR, CH4MMR,   &
     & C11MMR, C12MMR, O2MMR, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR,&
     & IS_NCOL, DT_BARE, DALB_BARE_WET, PEN_RAD_FRAC, SW_BETA

!- End of COMDECK declaration

      REAL :: CO2_MMR                ! CO2 concentration (if constant)

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
      ! Scaling factors to simulate inhomogeneous cloud.
      REAL, DIMENSION(NPD_CLOUD_COMPONENT) :: INHOM_CLOUD_SW
      REAL, DIMENSION(NPD_CLOUD_COMPONENT) :: INHOM_CLOUD_LW


! sza: add fraction standard diviation of inhomogeneous cloud
!      parameter for triple_clouds
!
      Real FW_STD

      ! Decorrelation pressure scale for large scale cloud
      REAL :: DP_CORR_STRAT
      ! Decorrelation pressure scale for convective cloud
      REAL :: DP_CORR_CONV

      COMMON  /RUN_Radiation/L_climat_aerosol, L_clim_aero_hgt,         &
     &  L_HadGEM1_Clim_Aero, FW_STD,                                    &
     &  A_SW_SEGMENTS,A_SW_SEG_SIZE,A_LW_SEGMENTS,A_LW_SEG_SIZE,CO2_MMR,&
     &  L_SEC_VAR,L_EqT,INHOM_CLOUD_SW,INHOM_CLOUD_LW,DP_CORR_STRAT,    &
     &  DP_CORR_CONV,AERO_BL_LEVELS, l_ovrlap, l_ccw_scav, l_cldtop_t_fix

      NAMELIST/RUN_Radiation/L_climat_aerosol, L_clim_aero_hgt,         &
     &  L_HadGEM1_Clim_Aero, FW_STD,                                    &
     &  A_SW_SEGMENTS,A_SW_SEG_SIZE,A_LW_SEGMENTS,A_LW_SEG_SIZE,CO2_MMR,&
     &  L_SEC_VAR,L_EqT,INHOM_CLOUD_SW,INHOM_CLOUD_LW,DP_CORR_STRAT,    &
     &  DP_CORR_CONV, AERO_BL_LEVELS,l_ovrlap,l_ccw_scav,l_cldtop_t_fix,&
      ! list of variables declared in <rad_com/rad_com.h> to be added:
     &  ALPHAM, ALPHAC, ALPHAB, DTICE,                                  &
     &  DT_BARE,DALB_BARE_WET,PEN_RAD_FRAC,SW_BETA,                     &
     &  N2OMMR, CH4MMR, C11MMR, C12MMR,                                 &
     &  O2MMR, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR, IS_NCOL

      ! Gravity wave drag:

      REAL :: KAY_GWAVE       ! Surface stress constant for GW drag

      REAL :: GWD_FRC         ! Critical Froude Number for 4A Scheme
       ! Maximum and minimum values for the STPH_RP scheme
       ! Gravity Wave drag
       REAL         :: GWD_FRC_max, GWD_FRC_min ! Max and min values
                       ! for the critical froud number
       REAL         :: KAY_GWAVE_max, KAY_GWAVE_min ! Max and min values
                       ! for the gravity wave constant

      REAL :: GWD_FSAT       ! Critical Froude number for wave breaking
                             ! used in the amplitude based saturation
                             ! test

      INTEGER :: SAT_SCHEME  ! Switch to determine whether to use
                             ! stress based    (sat_scheme=0) or
                             ! amplitude based (sat_scheme=1)
                             ! saturation test

      LOGICAL :: L_TAUS_SCALE  ! FROUDE NO. DEPENDENT SURFACE STRESS
      LOGICAL :: L_FIX_GWSATN  ! BUG FIXES TO GWSATN
      LOGICAL :: L_GWD_40km    ! Turn off orographic GWD above 40km
      LOGICAL :: L_USSP_OPAQUE ! OPAQUE LID FOR USSP SCHEME

      COMMON  /RUN_GWD/KAY_GWAVE, GWD_FRC, GWD_FSAT                     &
                      ,L_TAUS_SCALE, L_FIX_GWSATN, l_gwd_40km           &
                      ,L_USSP_OPAQUE, SAT_SCHEME                        &
                      ,GWD_FRC_max, GWD_FRC_min                         &
                      ,KAY_GWAVE_max, KAY_GWAVE_min

      NAMELIST/RUN_GWD/KAY_GWAVE, GWD_FRC, GWD_FSAT                     &
                      ,L_TAUS_SCALE, L_FIX_GWSATN, l_gwd_40km           &
                      ,L_USSP_OPAQUE, SAT_SCHEME                        &
                      ,GWD_FRC_max, GWD_FRC_min                         &
                      ,KAY_GWAVE_max, KAY_GWAVE_min

!------------------   End of Physics   ---------------------------------

      ! River Routing
      REAL :: RIVER_VEL       ! River velocity
      REAL :: RIVER_MCOEF     ! Meandering coefficient

      COMMON /RUN_RIVERS/ RIVER_VEL, RIVER_MCOEF
      NAMELIST /RUN_RIVERS/ RIVER_VEL, RIVER_MCOEF


      ! Aerosol Modelling

      INTEGER :: SO2_HIGH_LEVEL   ! Model level for chimney SO2 emiss
      INTEGER :: SOOT_HIGH_LEVEL  ! Model level for chimney soot emiss
      INTEGER :: BMASS_HIGH_LEVEL_1  ! Lowest and highest model levels
      INTEGER :: BMASS_HIGH_LEVEL_2  ! with biomass emissions.
      INTEGER :: OCFF_HIGH_LEVEL  ! Model level for chimney OCFF emiss

      ! improved nonhydrostatic  weights in tracer1 calculation
      LOGICAL :: L_tracer1_non_hydro

      COMMON  /RUN_Aerosol/                                             &
     & SO2_HIGH_LEVEL, SOOT_HIGH_LEVEL, BMASS_HIGH_LEVEL_1,             &
     & BMASS_HIGH_LEVEL_2, OCFF_HIGH_LEVEL, L_tracer1_non_hydro
      NAMELIST/RUN_Aerosol/                                             &
     & SO2_HIGH_LEVEL, SOOT_HIGH_LEVEL, BMASS_HIGH_LEVEL_1,             &
     & BMASS_HIGH_LEVEL_2, OCFF_HIGH_LEVEL, L_tracer1_non_hydro

!     Declarations for UKCA sub-model

      LOGICAL :: L_ukca_chem      ! True when UKCA chemistry is on
      LOGICAL :: L_ukca_family    ! True when using family chemistry
      LOGICAL :: L_ukca_advh2o    ! True when advecting H2O tracer
      LOGICAL :: L_ukca_phot2d    ! True when using 2D photolysis
      LOGICAL :: L_ukca_fastj     ! True when using Fastj photolysis
      LOGICAL :: L_ukca_trop      ! True for tropospheric chemistry
      LOGICAL :: L_ukca_tropisop  ! True for trop chemistry + isoprene
      LOGICAL :: L_ukca_strat     ! True for strat+reduced trop chemistry
      LOGICAL :: L_ukca_strattrop ! True for std strat+trop chemistry
      LOGICAL :: L_ukca_stratcfc  ! True for extended strat chemistry
      LOGICAL :: L_ukca_wachem    ! True for whole atm chemistry
      LOGICAL :: L_ukca_aerchem   ! True for trop+aerosol chemistry
      LOGICAL :: L_ukca_user      ! True for user-defined chemistry
      LOGICAL :: L_ukca_mode      ! True for UKCA-MODE aerosol scheme
      LOGICAL :: L_ukca_dust      ! True for UKCA dust scheme
      LOGICAL :: L_ukca_rnpb      ! True for Radon/Lead scheme
      LOGICAL :: L_ukca_o3budget  ! True for UKCA o3 budget
      LOGICAL :: L_ukca_qch4inter ! True for interact wetland CH4 ems
      LOGICAL :: L_ukca_isopinter ! True for interact isoprene ems
      LOGICAL :: L_ukca_terpinter ! True for interact terpene ems
      LOGICAL :: L_ukca_budget2   ! True for UKCA budget2 calculations
      LOGICAL :: L_ukca_qf11f12mbr  ! True for CFC emissions
      LOGICAL :: L_ukca_useumuivals ! True when using UMUI CFC values
      LOGICAL :: L_ukca_nat_sedi  ! True for NAT sedimentation
      LOGICAL :: L_ukca_het_psc   ! True for Het/PSC chemistry
      LOGICAL :: L_ukca_h2o_feedback ! True for H2O feedback from chem
      LOGICAL :: L_ukca_clbrcons  ! True for Cl/Br conservation
      LOGICAL :: L_ukca_usero3    ! True when using RO3 in photolysis
      LOGICAL :: L_ukca_useco3    ! True when using clim O3 in photolysis
      LOGICAL :: L_ukca_userelaxo3 ! True when using clim O3
!                                    above 0.3hPa in photolysis
      LOGICAL :: L_ukca_rado3   ! T when using UKCA O3 in radiation
      LOGICAL :: L_ukca_radch4  ! T when using UKCA CH4 in radiation
      LOGICAL :: L_ukca_radn2o  ! T when using UKCA N2O in radiation
      LOGICAL :: L_ukca_radf11  ! T when using UKCA CFC-11 in radn
      LOGICAL :: L_ukca_radf12  ! T when using UKCA CFC-12 in radn
      LOGICAL :: L_ukca_radf113 ! T when using UKCA CFC-113 in radn
      LOGICAL :: L_ukca_radf22  ! T when using UKCA HCFC-22 in radn
      LOGICAL :: L_ukca_radch2o ! T when using UKCA adv H2O in radn
      LOGICAL :: L_ukca_useoxid ! T when using UKCA oxidants in aer chem
      LOGICAL :: L_ukca_intdd   ! T when using interact dry deposition

! Tracers and chemistry integers:
      Integer :: jpctr          ! No. of chemical tracers
      Integer :: jpspec         ! No. of chemical species
      Integer :: jpbk           ! No. of bimolecular reactions
      Integer :: jptk           ! No. of termolecular reactions
      Integer :: jppj           ! No. of photolytic reactions
      Integer :: jphk           ! No. of heterogonous reactions
      Integer :: jpnr           ! jpbk + jptk + jppj + jphk
      Integer :: jpdd           ! No. of dry deposited species
      Integer :: jpdw           ! No. of wet deposited species
      Integer :: jpeq           ! No. of species in aqueous phase

! UKCA_MODE control features:
      LOGICAL :: L_ukca_primsu   ! T for primary sulphate aerosol emissions
      LOGICAL :: L_ukca_primss   ! T for primary sea-salt aerosol emissions
      LOGICAL :: L_ukca_primbcoc ! T for primary BC and OC aerosol emissions
      LOGICAL :: L_ukca_sedi     ! T for aerosol sedimentation
      LOGICAL :: L_ukca_nucl     ! T for aerosol sedimentation

      INTEGER :: i_mode_setup    ! Defines MODE aerosol scheme
      INTEGER :: i_mode_sizeprim ! Defines MODE size parameters
      INTEGER :: i_mode_nucscav  ! Defines nucleation scavenging method
      INTEGER :: i_mode_ddepaer  ! Defines dry deposition method
      INTEGER :: i_mode_nzts     ! No. of substeps for nucleation/sedimentation

      REAL :: mode_actdryr       ! Activation dry radius for MODE

      COMMON  /RUN_UKCA/ L_ukca_chem, L_ukca_family, L_ukca_advh2o,     &
     &         L_ukca_phot2d, L_ukca_fastj, L_ukca_trop,                &
     &         L_ukca_tropisop, L_ukca_strat, L_ukca_strattrop,         &
     &         L_ukca_stratcfc, L_ukca_wachem, L_ukca_aerchem,          &
     &         L_ukca_user, L_ukca_mode, L_ukca_dust,                   &
     &         L_ukca_rnpb, L_ukca_o3budget, L_ukca_qch4inter,          &
     &         L_ukca_isopinter, L_ukca_terpinter, L_ukca_budget2,      &
     &         L_ukca_qf11f12mbr, L_ukca_useumuivals, L_ukca_nat_sedi,  &
     &         L_ukca_het_psc, L_ukca_h2o_feedback, L_ukca_clbrcons,    &
     &         L_ukca_usero3, L_ukca_useco3, L_ukca_userelaxo3,         &
     &         L_ukca_rado3, L_ukca_radch4, L_ukca_radn2o,              &
     &         L_ukca_radf11, L_ukca_radf12, L_ukca_radf113,            &
     &         L_ukca_radf22, L_ukca_radch2o, L_ukca_useoxid,           &
     &         L_ukca_intdd, L_ukca_primsu, L_ukca_primss,              &
     &         L_ukca_primbcoc, L_ukca_sedi, L_ukca_nucl,               &
     &         jpctr, jpspec, jpbk, jptk, jppj, jphk, jpnr,             &
     &         jpdd, jpdw, jpeq,                                        &
     &         i_mode_setup, i_mode_sizeprim, i_mode_nucscav,           &
     &         i_mode_ddepaer, i_mode_nzts, mode_actdryr

      NAMELIST/RUN_UKCA/ L_ukca_chem, L_ukca_family, L_ukca_advh2o,     &
     &         L_ukca_phot2d, L_ukca_fastj, L_ukca_trop,                &
     &         L_ukca_tropisop, L_ukca_strat, L_ukca_strattrop,         &
     &         L_ukca_stratcfc, L_ukca_wachem, L_ukca_aerchem,          &
     &         L_ukca_user, L_ukca_mode, L_ukca_dust,                   &
     &         L_ukca_rnpb, L_ukca_o3budget, L_ukca_qch4inter,          &
     &         L_ukca_isopinter, L_ukca_terpinter, L_ukca_budget2,      &
     &         L_ukca_qf11f12mbr, L_ukca_useumuivals, L_ukca_nat_sedi,  &
     &         L_ukca_het_psc, L_ukca_h2o_feedback, L_ukca_clbrcons,    &
     &         L_ukca_usero3, L_ukca_useco3, L_ukca_userelaxo3,         &
     &         L_ukca_rado3, L_ukca_radch4, L_ukca_radn2o,              &
     &         L_ukca_radf11, L_ukca_radf12, L_ukca_radf113,            &
     &         L_ukca_radf22, L_ukca_radch2o, L_ukca_useoxid,           &
     &         L_ukca_intdd, L_ukca_primsu, L_ukca_primss,              &
     &         L_ukca_primbcoc, L_ukca_sedi, L_ukca_nucl,               &
     &         jpctr, jpspec, jpbk, jptk, jppj, jphk, jpnr,             &
     &         jpdd, jpdw, jpeq,                                        &
     &         i_mode_setup, i_mode_sizeprim, i_mode_nucscav,           &
     &         i_mode_ddepaer, i_mode_nzts, mode_actdryr

      INTEGER :: Instability_diagnostics  ! >0 if wanted, 0 otherwise
      INTEGER :: print_step    ! To control diagnostic printing interval
      INTEGER :: diag_interval ! diagnostic printing sampling frequency
      INTEGER :: norm_lev_start ! start level for norm diagnostics
      INTEGER :: norm_lev_end   ! end level for norm diagnostics
      INTEGER :: first_norm_print ! first timestep for norm printing
      INTEGER :: rpemax ! array size needed for diagnostic printing
      INTEGER :: rpemin ! array size needed for diagnostic printing
      INTEGER :: rpesum ! array size needed for diagnostic printing
      INTEGER :: ipesum ! array size needed for diagnostic printing
      INTEGER :: time_theta1_min      ! Timestep of min level 1 theta
      INTEGER :: time_w_max(max_model_levels) ! Timestep of max w
      INTEGER :: time_div_max(max_model_levels) ! Timestep of max div
      INTEGER :: time_div_min(max_model_levels) ! Timestep of min div
      INTEGER :: time_lapse_min(max_model_levels) ! Timestep of min
      INTEGER :: time_max_shear(max_model_levels) !Timestep max shear
      INTEGER :: time_max_wind(max_model_levels) ! Timestep of max wind
      INTEGER :: time_KE_max(max_model_levels) ! Timestep of max KE
      INTEGER :: time_KE_min(max_model_levels) ! Timestep of min KE
      INTEGER :: time_noise_max(max_model_levels) ! Timestep of max

      REAL:: frictional_timescale(max_model_levels) ! For idealised case
      REAL :: w_print_limit ! w Threshold for diagnostic printing
      REAL :: w_conv_limit  ! w Threshold for limiting convection
      REAL :: tropics_deg  ! define latitude for tropics
      REAL :: min_theta1_run                  ! Min theta level 1
      REAL :: dtheta1_run   ! Largest -ve delta theta at min theta1
      REAL :: max_w_run(max_model_levels)  ! Max w at a level
      REAL :: max_div_run(max_model_levels) ! Max divergence at a level
      REAL :: min_div_run(max_model_levels) ! Min divergence at a level
      REAL :: min_lapse_run(max_model_levels) ! Min dtheta/dz at a level
      REAL :: max_shear_run(max_model_levels) ! Max shear at a level
      REAL :: max_wind_run(max_model_levels) ! Max wind at a level
      REAL :: max_KE_run(max_model_levels)   ! Max KE at a level
      REAL :: min_KE_run(max_model_levels)   ! Min KE at a level
      REAL :: max_noise_run(max_model_levels) ! Max noise at a level

!     Problem_number not set here  ! Now controlled by namelist input
!     Instability_diagnostics      ! Now controlled by namelist input
!     frictional_timescale         ! Now intitialised in SETCONA

      LOGICAL :: L_idealised_data
      LOGICAL :: L_simple_friction
      LOGICAL :: L_print_w
      LOGICAL :: L_print_div
      LOGICAL :: L_diag_print_ops     ! diagnostic prints for ops
      LOGICAL :: L_print_pe     ! print diagnostics on all pe's if true
      LOGICAL :: L_print_shear        ! wind diagnostic prints
      LOGICAL :: L_print_max_wind     ! wind diagnostic prints
      LOGICAL :: L_diag_L2norms  ! l2norm diagnostic prints
      LOGICAL :: L_diag_L2helm   ! l2norm diagnostic prints from solver
      LOGICAL :: L_diag_wind     ! wind diagnostic prints
      LOGICAL :: L_diag_noise    ! w diagnostic prints
      LOGICAL :: L_flush6        ! if T then flush buffers on failure
      LOGICAL :: L_inviscid  ! Use prognostic inviscid theta at surface
      LOGICAL :: L_deep        ! Deep/shallow switch
      LOGICAL :: L_hydrostatic ! Hydrostatic option
      LOGICAL :: L_geostrophic ! Geostrophic balance at start
      LOGICAL :: L_solver      ! Solver active if true
      LOGICAL :: L_trap_errors ! Error trapping
      LOGICAL :: L_uv_zero     ! Set winds to zero
      LOGICAL :: L_diag_print  ! Print diagnostics
      LOGICAL :: L_print_lapse ! Print lapse_rate diagnostics
      LOGICAL :: L_print_wmax ! Print max w diagnostic
      LOGICAL :: L_print_theta1 ! Print level 1 theta diagnostic

!     L_idealised_data             ! Now controlled by namelist input

!------------------   Diagnostics:   --------------------------------

      COMMON  /RUN_Diagnostics/                                         &
     &  L_diag_print, L_diag_print_ops, L_flush6, print_step,           &
     &  diag_interval,w_print_limit, rpemax, rpemin, ipesum, rpesum,    &
     &  L_print_w, L_print_wmax,L_print_lapse, L_print_theta1,          &
     &  L_print_div, L_diag_wind, L_print_shear, L_print_max_wind,      &
     &  L_print_pe, L_diag_noise, L_diag_L2norms, L_diag_L2helm,        &
     &  norm_lev_start, norm_lev_end, first_norm_print,                 &
     &  max_w_run, min_theta1_run, dtheta1_run,                         &
     &  max_div_run, min_div_run, min_lapse_run,                        &
     &  max_shear_run, max_wind_run,                                    &
     &  max_noise_run, max_KE_run, min_KE_run,                          &
     &  time_KE_max, time_KE_min,                                       &
     &  time_w_max, time_div_max, time_div_min, time_lapse_min,         &
     &  time_max_shear, time_max_wind,                                  &
     &  time_theta1_min, time_noise_max

!------------------   Dynamics:   --------------------------------------

      ! Generalised dynamics:

      LOGICAL :: L_Physics   ! T: physics to be included
      LOGICAL :: L_Run_With_Physics2   ! T: physics2 to be included
      LOGICAL :: L_perturb_IC_theta    ! T: perturb theta on ts1   
      LOGICAL :: L_Backwards ! F: Integrate backwards without physics.
      LOGICAL :: L_Primitive ! T: Primitive Equation model
                   ! F: Semi-Geostrophic model (no dynamics)
      LOGICAL :: L_dry         ! T: run with no moisture
      LOGICAL :: L_adjust_wet  ! T: perform simple dry&moist adjustment
      LOGICAL :: L_mix_ratio

      LOGICAL :: L_LBC_balance  ! T: impose vertically balanced Exners
                                !    and rho in LBCs (set w=0)
                                ! F: leave Exners alone

      ! Switch for new (correct) lbcs developed in 2008
      LOGICAL :: L_lbc_new       ! .true. for new lbc switch

      LOGICAL :: L_free_slip    ! .true. for free-slip lower boundary

      ! T: reset polar values to mean value every polar_reset_timesteps
      LOGICAL, PARAMETER :: L_polar_reset = .false.

      ! interval for resetting polar to mean
      INTEGER, PARAMETER :: polar_reset_timesteps = 1

      INTEGER :: IntRand_seed

      ! Time weights for integration scheme:
      REAL :: alpha_1
      REAL :: alpha_2
      REAL :: alpha_3
      REAL :: alpha_4

!  variable resolution control
      LOGICAL :: L_regular  ! true if NOT variable resolution
      LOGICAL :: L_qwaterload  ! true if using adding waterloading terms

      INTEGER :: lam_var    ! number of variable res. lambda intervals
      INTEGER :: phi_var    ! number of variable res. phi intervals

      REAL :: var_ratio  ! grid-stretcing ratio for var grid
      REAL :: lam_ratio  ! scaling of original grid to high res grid
      REAL :: phi_ratio  ! scaling of original grid to high res grid
      REAL :: lam_frac  ! proportion of reg. points in West of domain
      REAL :: phi_frac  ! proportion of reg. points in South of domain
      REAL :: lambda_p_end
      REAL :: phi_p_end
      REAL :: lambda_u_end
      REAL :: phi_v_end
      REAL :: dlambda_p_end
      REAL :: dphi_p_end
      REAL :: dlambda_u_end
      REAL :: dphi_v_end

      ! GCR scheme= Generalized conjugate residual elliptic equation
      ! solver:
      LOGICAL :: GCR_use_tol_abs
      LOGICAL :: GCR_zero_init_guess
      LOGICAL :: GCR_use_residual_Tol

      ! T then use full equation on RHS on second and subsequent ADI
      ! timesteps
      LOGICAL :: GCR_adi_add_full_soln

      LOGICAL :: L_gcr_fast_x ! use faster non reproducible code
! Speed up solver convergence when dynamics-physics cycling
! is used.
       Logical :: L_GCR_cycle_opt

      LOGICAL :: L_interp_depart ! interpolate to find u,v departure pts

! ----------------------- include file: GCR_ITER_DIM -------------------
! Description: Introduce integer parameters to:
!              (1) declare the maximum allowed number of dynamics 
!                  cycles (timestep iterations)
!              (2) declare the number of GCR Iterations analysis steps
!
! Code Author: M. Diamantakis
! Current Code Owner: A. Malcom
!
      INTEGER, PARAMETER :: MAX_NUMCYCLES = 10 ! Max number of 
                                               ! dynamics cycles
      INTEGER, PARAMETER :: GCR_ANAL_STEPS = 3 ! Number of GCR Iterations
                                               ! analysis step

      INTEGER :: GCR_max_iterations
      INTEGER :: GCR_diagnostics    !RR replace by STASHflag?
      INTEGER :: GCR_its_avg_step(GCR_anal_steps)   
                                             ! GCR Iterations analysis 
                                             ! steps
      INTEGER, DIMENSION(max_numcycles) :: GCR_its_switch
                                             ! GCR Iterations analysis 
                                             ! switch
      INTEGER, DIMENSION(max_numcycles) :: GCR_max_its  
                                             ! GCR Max iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_min_its  
                                             ! GCR Min iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_max_time 
                                             ! Timestep number
                                             ! GCR Max iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_min_time 
                                             ! Timestep number
                                             ! GCR Min iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_sum_its  
                                             ! Sum iterations
                                             ! over test period
      INTEGER :: GCR_precon_option
      INTEGER :: GCR_n_ADI_pseudo_timesteps

      ! After how many iterations do we restart
      INTEGER :: GCR_Restart_value

      REAL :: GCR_tol_res
      REAL :: GCR_tol_abs
      REAL :: GCR_ADI_pseudo_timestep
      REAL :: G_term_tol
!      G_term_tol             ! Now controlled by namelist input
!
!  Parameters for dynamics-physics cycling
       INTEGER :: NumCycles
       LOGICAL :: L_new_tdisc
       REAL    :: extrp_weight
       REAL    :: GCR_tol_abs2
       REAL    :: GCR_tol_res2
       REAL    :: alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2
! Parameter for fully-interpolating theta option
       LOGICAL :: L_fint_theta
! Number of physics2 substeps (default is 1)
       INTEGER :: Num_Substeps
! Logical switch to enable phys2 substepping
       Logical :: L_phys2_substep
      COMMON  /RUN_Dyn/                                                 &
     &  L_idealised_data,                                               &
     &  L_Physics, L_phys2_substep, L_Backwards, L_Primitive,           &
     &  L_GCR_cycle_opt, L_Run_With_Physics2, L_perturb_IC_theta,       &
     &  L_mix_ratio, L_new_tdisc, L_fint_theta, IntRand_seed,           &
     &  L_free_slip,                                                    &
     &  L_simple_friction, L_inviscid, L_deep, L_hydrostatic,           &
     &  L_dry, L_adjust_wet, L_LBC_balance, L_lbc_new,                  &
     &  L_geostrophic, L_solver, L_trap_errors, L_uv_zero,              &
     &  L_regular, L_qwaterload, n_rims_to_do, lam_var, phi_var,        &
     &  var_ratio, lam_ratio, phi_ratio, lam_frac, phi_frac,            &
     &  lambda_p_end,  phi_p_end, lambda_u_end,  phi_v_end,             &
     &  dlambda_p_end,  dphi_p_end, dlambda_u_end,  dphi_v_end,         &
     &  alpha_1,alpha_2,alpha_3,alpha_4,                                &
     &  alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2,                     &
     &  Num_Substeps, NumCycles,                                        &
     &  extrp_weight,                                                   &
     &  GCR_diagnostics, GCR_max_its, GCR_min_its,                      &
     &  GCR_max_time, GCR_min_time,                                     &
     &  GCR_sum_its, GCR_its_switch, GCR_its_avg_step,                  &
     &  GCR_use_tol_abs,GCR_zero_init_guess,GCR_use_residual_Tol,       &
     &  GCR_adi_add_full_soln,L_gcr_fast_x,                             &
     &  L_interp_depart,                                                &
     &  GCR_max_iterations, GCR_precon_option,                          &
     &  GCR_n_ADI_pseudo_timesteps,GCR_Restart_value,GCR_tol_res,       &
     &  GCR_tol_res2,GCR_tol_abs,GCR_tol_abs2,GCR_ADI_pseudo_timestep,  &
     &  G_term_tol
      NAMELIST/RUN_Dyn/                                                 &
     &  L_Physics,L_phys2_substep,L_Backwards,L_Primitive,              &
     &  L_GCR_cycle_opt, L_Run_With_Physics2, L_perturb_IC_theta,       &
     &  L_mix_ratio, L_new_tdisc, L_fint_theta, IntRand_seed,           &
     &  L_free_slip, L_dry, L_adjust_wet,                               &
     &  L_LBC_balance, L_lbc_new,                                       &
     &  L_regular, L_qwaterload, n_rims_to_do, lam_var, phi_var,        &
     &  var_ratio, lam_ratio, phi_ratio, lam_frac, phi_frac,            &
     &  lambda_p_end,  phi_p_end, lambda_u_end,  phi_v_end,             &
     &  dlambda_p_end,  dphi_p_end, dlambda_u_end,  dphi_v_end,         &
     &  alpha_1,alpha_2,alpha_3,alpha_4,                                &
     &  alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2,                     &
     &  Num_Substeps, NumCycles,                                        &
     &  extrp_weight,                                                   &
     &  GCR_use_tol_abs,GCR_zero_init_guess,GCR_use_residual_Tol,       &
     &  GCR_adi_add_full_soln,L_gcr_fast_x,                             &
     &  L_interp_depart,                                                &
     &  GCR_max_iterations,GCR_diagnostics,GCR_precon_option,           &
     &  GCR_n_ADI_pseudo_timesteps,GCR_Restart_value,GCR_tol_res,       &
     &  GCR_tol_res2,GCR_tol_abs,GCR_tol_abs2,GCR_ADI_pseudo_timestep,  &

     &  G_term_tol, GCR_its_avg_step

! Semi-Lagrangian advection control options:


      INTEGER, PARAMETER :: Number_SL_choices = 4
      INTEGER, PARAMETER :: Theta_SL = 1
      INTEGER, PARAMETER :: moist_SL = 2
      INTEGER, PARAMETER :: Wind_SL = 3
      INTEGER, PARAMETER :: rho_SL = 4

      ! T if conservation to be enforced.
      LOGICAL:: L_conserv(Number_SL_choices)

      ! T if monotonicity to be enforced.
      LOGICAL:: L_mono(Number_SL_choices)

      ! T if high order interpolation scheme to be used.
      LOGICAL:: L_high(Number_SL_choices)

      ! T if high order scheme to be used in Ritchie routine.
      LOGICAL:: L_Ritchie_high

      ! T if monotone scheme to be used in Ritchie routine.
      LOGICAL:: L_Ritchie_mono

      ! T then only perform vector co-ordinate geometry in 2d.
      LOGICAL:: L_2d_sl_geometry


      ! T: sl code bit repoducible with any sensible halo size
      LOGICAL:: L_sl_halo_reprod

      ! improved nonhydrostatic conservation  of q, qcl and qcf
      LOGICAL:: L_moist_nonhydro_conserve

      ! Run tracer advection code with conservation on          
      LOGICAL :: L_conserve_tracers

      ! a code saying which high monotone scheme to use.
      ! 1 = tensor tri-cubic lagrange (j,i,k)
      ! no other options available at present.
      INTEGER :: high_order_scheme(Number_SL_choices)

      ! a code saying which  monotone scheme to use.
      ! 1 = tri-linear order (j,i,k)
      ! no other options available at present.
      INTEGER :: monotone_scheme(Number_SL_choices)
      REAL :: thmono_height       ! top for monotone fully-interp theta
      INTEGER :: thmono_levels ! levels for monotone fully-interp theta
!                 increments used to limit non-interp increments

      ! a code saying which high  order scheme to use in Ritchie routine
      INTEGER :: Ritchie_high_order_scheme

      ! a code saying which monotone scheme to use in Ritchie routine.
      INTEGER :: Ritchie_monotone_scheme

      INTEGER :: Depart_scheme     ! which departure point scheme to use

      ! for the chosen departure point scheme how many iterations/terms
      ! to use
      INTEGER :: Depart_order

      ! used in interpolation code:
      ! number of levels either side of default level to search.
      INTEGER :: interp_vertical_search_tol

      ! Set of look-up tables for searching on variable grids
      INTEGER :: look_lam(max_look)  ! lambda p search
      INTEGER :: look_phi(max_look)  ! phi p search

      ! Set of minumum search grid lengths for variable grids
      REAL ::  recip_dlam   ! smallest delta_lambda p
      REAL ::  recip_dphi   ! smallest delta_phi p

      INTEGER :: halo_lam  ! halo for lamp look-up table
      INTEGER :: halo_phi  ! halo for phip look-up table


      COMMON  /RUN_SL/L_conserv,L_mono,L_high,L_Ritchie_high,           &
     &  L_Ritchie_mono, thmono_height, thmono_levels, L_2d_sl_geometry, &
     &  L_sl_halo_reprod,high_order_scheme,monotone_scheme,             &
     &  Ritchie_high_order_scheme,Ritchie_monotone_scheme,              &
     &  Depart_scheme,Depart_order,interp_vertical_search_tol,          &
     &  recip_dlam, recip_dphi, halo_lam, halo_phi,                     &
     &  look_lam, look_phi,                                             &
     &  L_moist_nonhydro_conserve, L_conserve_tracers,                  &
     &  Instability_diagnostics

      NAMELIST/RUN_SL/L_conserv,L_mono,L_high,L_Ritchie_high,           &
     &  L_Ritchie_mono, thmono_height, L_2d_sl_geometry,                &
     &  L_sl_halo_reprod,high_order_scheme,monotone_scheme,             &
     &  Ritchie_high_order_scheme,Ritchie_monotone_scheme,              &
     &  Depart_scheme,Depart_order,interp_vertical_search_tol,          &
     &  L_moist_nonhydro_conserve, L_conserve_tracers

      ! Diffusion control options:
      LOGICAL :: L_diffusion
      LOGICAL :: L_upper_ramp   ! ramp upper-level diffusion
      LOGICAL :: L_adjust_theta ! activate convective adjustment
      LOGICAL :: L_vdiff_uv    ! activate targeted diffusion uv
      LOGICAL :: L_filter    ! activate polar filter or diffusion
      LOGICAL :: L_filter_incs    ! activate polar filter of incs
      LOGICAL :: L_pftheta   ! activate polar filter of theta
      LOGICAL :: L_pfuv      ! activate polar filter of u,v
      LOGICAL :: L_pfw       ! activate polar filter of w
      LOGICAL :: L_pofil_new  ! activate new polar filter
      LOGICAL :: L_pfexner  ! activate polar filter of Exner pressure
      LOGICAL :: L_pofil_hadgem2 ! run with HadGEM2 polar filter setting   
      LOGICAL :: L_diff_exner  ! activate diffusion of Exner pressure
      LOGICAL :: L_pfcomb   ! combined polar filter/diffusion active
      LOGICAL :: L_sponge   ! activate lateral boundariessponge zones
      LOGICAL :: L_pfincs    ! activate polar filter of incs
      LOGICAL :: L_diff_thermo ! horiz. diffusion of theta
      LOGICAL :: L_diff_wind   ! horiz. diffusion of u, v
      LOGICAL :: L_diff_w      ! horiz. diffusion of w
      LOGICAL :: L_diff_incs   ! horiz. diffusion of increments
      LOGICAL :: L_diff_auto   ! UM calculates diffusion parameters
      LOGICAL :: L_tardiff_q     ! activate targeted diffusion q
      LOGICAL :: L_diff_ctl      ! general diffusion control
! L_diff_ctl == L_diffusion .or. L_cdiffusion .or. L_vertical_diffusion
!               .or. L_divdamp .or.  L_tardiff_q .or. L_diag_print
      LOGICAL :: L_cdiffusion
      LOGICAL :: L_subfilter_horiz  ! subgrid turbulence in horizontal
      LOGICAL :: L_subfilter_vert   ! subgrid turbulence in vertical  
      LOGICAL :: L_subfilter_blend  ! blend diffusion coeffs

      ! value * del^2: order=2
      INTEGER  :: diffusion_order_thermo(max_model_levels)

      ! gives del^4 diffusion
      INTEGER  :: diffusion_order_wind(max_model_levels)
      INTEGER  :: diffusion_order_w(max_model_levels)
      INTEGER  :: diffusion_order_q(max_model_levels)
      INTEGER :: u_begin(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: u_end(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_begin(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_end(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: u_sweeps(max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_sweeps(max_121_rows) ! Sweep control on 121 filter
      INTEGER :: max_sweeps ! Max sweeps wanted for 121 filter
      INTEGER :: global_u_filter ! Sweep control on 121 filter
      INTEGER :: global_v_filter ! Sweep control on 121 filter
      INTEGER :: diff_order_thermo   ! diffusion order for theta
      INTEGER :: diff_order_wind ! diffusion order for winds
      INTEGER :: diff_timescale_thermo  ! diffusion timescale for theta
      INTEGER :: diff_timescale_wind    ! diffusion timescale for wind
      INTEGER :: vdiffuv_timescale ! diffusion e-folding timesteps
      INTEGER :: vdiffuv_start ! start level  targeted diffusion
      INTEGER :: vdiffuv_end ! end level targeted diffusion
      INTEGER :: adjust_theta_start ! start level convective adjustment
      INTEGER :: adjust_theta_end   ! end levelconvective adjustment
      INTEGER :: top_filt_start ! start level upper-level diffusion
      INTEGER :: top_filt_end   ! end level upper-level diffusion
      INTEGER :: sponge_ew   ! left/right boundaries sponge zone width
      INTEGER :: sponge_ns   ! north/south boundaries sponge zone width
      INTEGER :: sponge_power ! sponge zone weighting order
      INTEGER :: tardiffq_test ! test level test w targetted diffusion
      INTEGER :: tardiffq_start ! start level test w targetted diffusion
      INTEGER :: tardiffq_end ! end level test w targetted diffusion
      INTEGER :: turb_startlev_horiz   ! 1st lev for horiz subgrid turb
      INTEGER :: turb_endlev_horiz    ! last lev for horiz subgrid turb
      INTEGER :: turb_startlev_vert   ! 1st lev for vert subgrid turb
      INTEGER :: turb_endlev_vert     ! last lev for vert subgrid turb

      !  level - assume surfaces are horizontal
      INTEGER  :: horizontal_level
      INTEGER  ::  tar_horizontal ! steep slope test targeted diffusion

      REAL :: diffusion_coefficient_thermo(max_model_levels)
      REAL :: diffusion_coefficient_wind(max_model_levels)
      REAL :: diffusion_coefficient_w(max_model_levels)
      REAL :: diffusion_coefficient_q(max_model_levels)
      REAL :: vdiffuv_test ! test to activate shear diffusion of u,v
      REAL :: vdiffuv_factor ! vertical diffusion coeff for shear diff
      REAL :: scale_ratio ! Pass control on 121 filter
      REAL :: ref_lat_deg ! Reference lat for auto diffusion
      REAL :: top_diff    !  upper-level diffusion coefficient
      REAL :: up_diff_scale    !  upper-level diffusion ramping factor

      REAL :: up_diff(max_updiff_levels) ! upper-level diffusion coeff
      REAL :: sponge_wts_ew(max_sponge_width) ! sponge weights
      REAL :: sponge_wts_ns(max_sponge_width) ! sponge weights

      REAL :: diff_coeff_ref ! EW diffusion coefficient at polar cap
      REAL :: diff_coeff_thermo  ! NS theta diffusion coeff
      REAL :: diff_coeff_wind    ! NS u,v diffusion coeff
      REAL :: diff_coeff_phi    ! North-South diffusion coefficient
!   reference latitudes for filtering and diffusion
      REAL :: polar_cap  ! Apply 1-2-1 filter polewards
      REAL :: tardiffq_factor ! targeted diffusion coefficient
      REAL :: diff_factor
      REAL :: mix_factor

      ! Divergence damping control options:
      LOGICAL :: L_divdamp

      REAL :: div_damp_coefficient(max_model_levels)

      ! Polar filter control options:
      LOGICAL :: L_polar_filter
      ! T: use polar filter to filter increment
      LOGICAL :: L_polar_filter_incs

      REAL :: polar_filter_north_lat_limit
      REAL :: polar_filter_south_lat_limit
      REAL :: polar_filter_coefficient

      ! amount in radians to increment start latitude by per sweep
      REAL :: polar_filter_step_per_sweep

      ! max latitude at which filter can star
      REAL :: polar_filter_lat_limit

      ! number of sweeps of filter to do
      INTEGER :: polar_filter_n_sweeps

      ! Vertical Diffusion control options:
      LOGICAL :: L_vertical_diffusion
      LOGICAL :: L_ramp

      INTEGER :: level_start_wind
      INTEGER :: level_stop_wind
      INTEGER :: level_start_theta
      INTEGER :: level_stop_theta
      INTEGER :: level_start_q
      INTEGER :: level_stop_q

      REAL :: vert_diffusion_coeff_wind
      REAL :: vert_diffusion_coeff_theta
      REAL :: vert_diffusion_coeff_q
      REAL :: ramp_lat_radians

      ! Moisture resetting control options:
      LOGICAL :: l_qpos         ! logical to run qpos code

      ! true to do Method 2 false to do Method 1
      LOGICAL :: l_q_pos_local


      REAL :: qlimit         !    lowest allowed value of q

      COMMON  /RUN_Diffusion/L_diffusion,diffusion_order_thermo,        &
     &  diffusion_order_wind,diffusion_order_w,diffusion_order_q,       &
     &  diffusion_coefficient_thermo,diffusion_coefficient_wind,        &
     &  diffusion_coefficient_w,diffusion_coefficient_q,                &
     &  diff_order_thermo, diff_timescale_thermo, diff_coeff_thermo,    &
     &  diff_order_wind, diff_timescale_wind, diff_coeff_wind,          &
     &  horizontal_level, tar_horizontal, tropics_deg,                  &
     &  L_cdiffusion,L_ramp, ramp_lat_radians,                          &
     &  L_divdamp,div_damp_coefficient,                                 &
     &  L_polar_filter,L_polar_filter_incs,                             &
     &  polar_filter_north_lat_limit,polar_filter_south_lat_limit,      &
     &  polar_filter_coefficient,polar_filter_step_per_sweep,           &
     &  polar_filter_lat_limit,polar_filter_n_sweeps,                   &
     &  l_qpos,  l_q_pos_local, qlimit,                                 &
     &  L_diff_ctl, L_tardiff_q, w_conv_limit,                          &
     &  tardiffq_factor, tardiffq_test, tardiffq_start, tardiffq_end,   &
     &  L_subfilter_horiz, L_subfilter_vert, L_subfilter_blend,         &
     &  diff_factor, mix_factor, turb_startlev_horiz, turb_endlev_horiz,&
     &  turb_startlev_vert, turb_endlev_vert,                           &
     &  L_vertical_diffusion, L_pofil_hadgem2,                          &
     &  L_pofil_new, L_pfcomb, L_pftheta, L_pfuv, L_pfw, L_pfincs,      &
     &  L_pfexner, L_filter, L_filter_incs, L_diff_exner, L_diff_auto,  &
     &  L_diff_incs, L_diff_thermo, L_diff_wind, L_diff_w,              &
     &  L_vdiff_uv, L_adjust_theta,                                     &
     &  level_start_wind, level_stop_wind,                              &
     &  level_start_theta, level_stop_theta,                            &
     &  level_start_q, level_stop_q,                                    &
     &  polar_cap, scale_ratio, ref_lat_deg, max_sweeps,                &
     &  L_upper_ramp, top_filt_start, top_filt_end, top_diff,           &
     &  up_diff, up_diff_scale, sponge_wts_ew, sponge_wts_ns,           &
     &  u_begin, u_end, v_begin, v_end, u_sweeps, v_sweeps,             &
     &  global_u_filter, global_v_filter,                               &
     &  diff_coeff_ref, diff_coeff_phi, vdiffuv_timescale,              &
     &  vdiffuv_test, vdiffuv_factor, vdiffuv_start, vdiffuv_end,       &
     &  adjust_theta_start, adjust_theta_end,                           &
     &  L_sponge, sponge_power, sponge_ew, sponge_ns,                   &
     &  vert_diffusion_coeff_wind, vert_diffusion_coeff_theta,          &
     &  vert_diffusion_coeff_q

      NAMELIST/RUN_Diffusion/L_diffusion,diffusion_order_thermo,        &
     &  diffusion_order_wind,diffusion_order_w,diffusion_order_q,       &
     &  diffusion_coefficient_thermo,diffusion_coefficient_wind,        &
     &  diffusion_coefficient_w,diffusion_coefficient_q,                &
     &  diff_order_thermo, diff_timescale_thermo,                       &
     &  diff_order_wind, diff_timescale_wind,                           &
     &  horizontal_level, tar_horizontal,                               &
     &  L_cdiffusion,L_ramp, ramp_lat_radians,                          &
     &  L_divdamp,div_damp_coefficient,                                 &
     &  L_polar_filter,L_polar_filter_incs,                             &
     &  polar_filter_north_lat_limit,polar_filter_south_lat_limit,      &
     &  polar_filter_coefficient,polar_filter_step_per_sweep,           &
     &  polar_filter_lat_limit,polar_filter_n_sweeps,                   &
     &  l_qpos,  l_q_pos_local, qlimit,                                 &
     &  L_diff_ctl, L_tardiff_q, w_conv_limit,                          &
     &  tardiffq_factor, tardiffq_test, tardiffq_start, tardiffq_end,   &
     &  L_subfilter_horiz, L_subfilter_vert, L_subfilter_blend,         &
     &  diff_factor, mix_factor, turb_startlev_horiz, turb_endlev_horiz,&
     &  turb_startlev_vert, turb_endlev_vert,                           &
     &  L_vertical_diffusion,                                           &
     &  L_pftheta, L_pfuv, L_pfw, L_pfincs, L_pfexner, L_pofil_hadgem2, &
     &  L_diff_incs, L_diff_thermo, L_diff_wind, L_diff_w,              &
     &  level_start_wind, level_stop_wind,                              &
     &  level_start_theta, level_stop_theta,                            &
     &  level_start_q, level_stop_q,                                    &
     &  L_pofil_new, L_diff_auto,                                       &
     &  diff_coeff_ref, polar_cap, scale_ratio, ref_lat_deg,            &
     &  max_sweeps, L_upper_ramp, up_diff_scale, top_diff,              &
     &  top_filt_start, top_filt_end, L_vdiff_uv, vdiffuv_timescale,    &
     &  vdiffuv_test, vdiffuv_factor, vdiffuv_start, vdiffuv_end,       &
     &  L_adjust_theta, adjust_theta_start, adjust_theta_end,           &
     &  L_sponge, sponge_power, sponge_ew, sponge_ns,                   &
     &  vert_diffusion_coeff_wind, vert_diffusion_coeff_theta,          &
     &  vert_diffusion_coeff_q,                                         &
     &  L_diag_print, L_diag_print_ops, L_print_pe,                     &
     &  L_print_w, L_print_wmax, L_print_lapse, L_print_theta1,         &
     &  L_print_div, L_diag_wind, L_print_shear, L_print_max_wind,      &
     &  L_diag_noise, L_diag_L2norms, L_diag_L2helm,                    &
     &  norm_lev_start, norm_lev_end, first_norm_print,                 &
     &  print_step, diag_interval, w_print_limit, L_Flush6

!------------------   Dynamical core   -------------------------------
! Suarez-Held variables
      REAL :: SuHe_newtonian_timescale_ka
      REAL :: SuHe_newtonian_timescale_ks
      REAL :: SuHe_pole_equ_deltaT
      REAL :: SuHe_static_stab
      REAL :: base_frictional_timescale
      REAL :: SuHe_sigma_cutoff
      REAL :: SuHe_level_weight(max_model_levels)
      REAL :: friction_level(max_model_levels)

      INTEGER :: SuHe_relax
      INTEGER :: SuHe_fric

      LOGICAL :: L_SH_Williamson

      COMMON/Run_Dyncore/                                               &
     &  SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,       &
     &  SuHe_pole_equ_deltaT, SuHe_static_stab,                         &
     &  base_frictional_timescale, SuHe_sigma_cutoff,                   &
     &  L_SH_Williamson, SuHe_relax, SuHe_fric,                         &
     &  SuHe_level_weight, frictional_timescale, friction_level

!------------------  Idealised model   ----------------------------

      INTEGER,PARAMETER:: max_num_profile_data = 100
      INTEGER,PARAMETER:: max_num_force_times = 100
      INTEGER,PARAMETER:: idl_max_num_bubbles = 3

! Idealised  variables
      REAL :: h_o
      REAL :: h_o_actual  ! height of growing mountain
      REAL :: h_o_per_step ! height change per step of growing mountain
      REAL :: lambda_fraction
      REAL :: phi_fraction
      REAL :: half_width_x
      REAL :: half_width_y
      REAL :: Witch_power
      REAL :: plat_size_x
      REAL :: plat_size_y
      REAL :: height_domain
      REAL :: delta_x
      REAL :: delta_y
      REAL :: big_factor
      REAL :: mag
      REAL :: first_theta_height
      REAL :: thin_theta_height
      REAL :: p_surface
      REAL :: theta_surface
      REAL :: dtheta_dz1(3)
      REAL :: height_dz1(3)
      REAL :: Brunt_Vaisala
      REAL :: u_in(4)
      REAL :: v_in(4)
      REAL :: height_u_in(3)
      REAL :: u_ramp_start
      REAL :: u_ramp_end
      REAL :: ujet_lat
      REAL :: ujet_width
      REAL :: t_horizfn_data(10)
      REAL :: q1
      REAL :: theta_ref(max_model_levels)
      REAL :: rho_ref(max_model_levels)
      REAL :: exner_ref(max_model_levels + 1)
      REAL :: q_ref(max_model_levels)
      REAL :: u_ref(max_model_levels)
      REAL :: v_ref(max_model_levels)
      REAL :: z_orog_print(0:max_model_levels)
      REAL :: f_plane
      REAL :: ff_plane
      REAL :: r_plane
      REAL :: zprofile_data(max_num_profile_data)
      REAL :: tprofile_data(max_num_profile_data)
      REAL :: qprofile_data(max_num_profile_data)
      REAL :: z_uvprofile_data(max_num_profile_data)
      REAL :: uprofile_data(max_num_profile_data)
      REAL :: vprofile_data(max_num_profile_data)
      REAL :: tforce_time_interval
      REAL :: qforce_time_interval
      REAL :: uvforce_time_interval
      REAL :: newtonian_timescale
      REAL :: z_tforce_data(max_num_profile_data)
      REAL :: z_qforce_data(max_num_profile_data)
      REAL :: z_uvforce_data(max_num_profile_data)
      REAL :: tforce_data(max_num_profile_data, max_num_force_times)
      REAL :: qforce_data(max_num_profile_data, max_num_force_times)
      REAL :: uforce_data(max_num_profile_data, max_num_force_times)
      REAL :: vforce_data(max_num_profile_data, max_num_force_times)
      REAL :: tforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: qforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: uforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: vforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: pforce_time_interval
      REAL :: p_surface_data(max_num_force_times)
      REAL :: perturb_factor
      REAL :: perturb_magnitude_t
      REAL :: perturb_magnitude_q
      REAL :: perturb_height(2)
      REAL :: orog_hgt_lbc
      REAL :: zprofile_orog
      REAL :: hf
      REAL :: cool_rate
      REAL :: IdlSurfFluxSeaParams(10) ! Idealised surface flux params
      REAL :: roughlen_z0m   
      REAL :: roughlen_z0h
      ! Idealised bubbles
      REAL :: idl_bubble_max(idl_max_num_bubbles) ! Bubble max amplitude
      REAL :: idl_bubble_height(idl_max_num_bubbles)  ! Bubble height
      REAL :: idl_bubble_width(idl_max_num_bubbles)   ! Bubble width
      REAL :: idl_bubble_depth(idl_max_num_bubbles)   ! Bubble depth
      ! Bubble x-offset, y-offset in normalised units (0:1)
      ! (0.5=domain centre)
      REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
      REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
      REAL :: DMPTIM, HDMP, ZDMP   ! Damping layer values
      REAL :: u_geo, v_geo         ! Geostrophic wind


      INTEGER :: n_rims_to_do   ! rim size for LAM
      INTEGER :: surface_type
      INTEGER :: grow_steps
      INTEGER :: grid_number
      INTEGER :: grid_flat
      INTEGER :: tprofile_number
      INTEGER :: qprofile_number
      INTEGER :: uvprofile_number
      INTEGER :: num_profile_data
      INTEGER :: num_uvprofile_data
      INTEGER :: t_horizfn_number
      INTEGER :: uv_horizfn_number
      INTEGER :: pforce_option
      INTEGER :: num_pforce_times
      INTEGER :: tforce_option
      INTEGER :: qforce_option
      INTEGER :: uvforce_option
      INTEGER :: num_tforce_levels
      INTEGER :: num_tforce_times
      INTEGER :: num_qforce_levels
      INTEGER :: num_qforce_times
      INTEGER :: num_uvforce_levels
      INTEGER :: num_uvforce_times
      INTEGER :: IdlSurfFluxSeaOption  ! Idealised surface flux option
      INTEGER :: first_constant_r_rho_level_new
      INTEGER :: big_layers
      INTEGER :: transit_layers
      INTEGER :: mod_layers
      INTEGER :: idl_bubble_option(idl_max_num_bubbles) ! Bubble option
      INTEGER :: idl_interp_option  ! Profile interpolation option
      INTEGER :: perturb_type

      LOGICAL :: L_initialise_data
      LOGICAL :: L_constant_dz
      LOGICAL :: L_trivial_trigs !.false. for Cartesian coords (lat=0.0)
      LOGICAL :: L_idl_bubble_saturate(idl_max_num_bubbles)
      LOGICAL :: L_fixed_lbcs
      LOGICAL :: L_fix_orog_hgt_lbc
      LOGICAL :: L_pressure_balance
      LOGICAL :: L_wind_balance
      LOGICAL :: L_rotate_winds
      LOGICAL :: L_polar_wind_zero
      LOGICAL :: L_vert_Coriolis
      LOGICAL :: L_rotating     ! .true. for Earth's rotation
      LOGICAL :: L_perturb      ! add random perturb. to surface theta
      LOGICAL :: L_code_test    ! User switch for testing code
      LOGICAL :: L_pforce
      LOGICAL :: L_baroclinic
      LOGICAL :: L_cyclone
      LOGICAL :: L_force
      LOGICAL :: L_force_lbc
      LOGICAL :: L_perturb_t
      LOGICAL :: L_perturb_q
      LOGICAL :: L_perturb_correlate_tq
      LOGICAL :: L_perturb_correlate_vert
      LOGICAL :: L_perturb_correlate_time
      LOGICAL :: L_damp      ! Logical for damping layer
      LOGICAL :: L_geo_for ! Logical for geostrophic wind forcing
      LOGICAL :: L_bomex     ! Logical for BOMEX set up
      LOGICAL :: L_spec_z0   ! specification of roughness length    

      COMMON  /RUN_Ideal/                                               &
     &  h_o, h_o_actual, h_o_per_step,                                  &
     &  lambda_fraction, phi_fraction, half_width_x, half_width_y,      &
     &  Witch_power, plat_size_x, plat_size_y,                          &
     &  height_domain, delta_x, delta_y, big_factor, mag,               &
     &  first_theta_height, thin_theta_height, p_surface,               &
     &  theta_surface, dtheta_dz1, height_dz1, Brunt_Vaisala,           &
     &  u_in, v_in, height_u_in, u_ramp_start, u_ramp_end, q1,          &
     &  ujet_lat, ujet_width,                                           &
     &  t_horizfn_number, t_horizfn_data, uv_horizfn_number,            &
     &  u_ref, v_ref, theta_ref, exner_ref, rho_ref, q_ref,             &
     &  z_orog_print, grow_steps,                                       &
     &  surface_type, grid_number, grid_flat,                           &
     &  tprofile_number, qprofile_number, uvprofile_number,             &
     &  num_profile_data, num_uvprofile_data,                           &
     &  tforce_option, qforce_option, uvforce_option,                   &
     &  num_tforce_levels, num_tforce_times,                            &
     &  num_qforce_levels, num_qforce_times,                            &
     &  num_uvforce_levels, num_uvforce_times,                          &
     &  L_pforce, pforce_option, num_pforce_times,                      &
     &  first_constant_r_rho_level_new,                                 &
     &  big_layers, transit_layers, mod_layers,                         &
     &  zprofile_data, tprofile_data, qprofile_data,                    &
     &  z_uvprofile_data, uprofile_data, vprofile_data,                 &
     &  tforce_time_interval, qforce_time_interval,                     &
     &  uvforce_time_interval, newtonian_timescale,                     &
     &  z_tforce_data, z_qforce_data, z_uvforce_data,                   &
     &  tforce_data, qforce_data, uforce_data, vforce_data,             &
     &  tforce_data_modlev, qforce_data_modlev,                         &
     &  uforce_data_modlev, vforce_data_modlev,                         &
     &  pforce_time_interval, p_surface_data,                           &
     &  L_initialise_data,                                              &
     &  L_perturb_t, perturb_magnitude_t,                               &
     &  L_perturb_q, perturb_magnitude_q,                               &
     &  L_perturb_correlate_tq,                                         &
     &  L_perturb_correlate_vert,                                       &
     &  L_perturb_correlate_time,                                       &
     &  perturb_type, perturb_height,                                   &
     &  L_constant_dz, L_polar_wind_zero,                               &
     &  L_wind_balance, L_rotate_winds,                                 &
     &  IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                     &
     &  L_spec_z0, roughlen_z0m, roughlen_z0h,                          &
     &  L_pressure_balance, L_vert_Coriolis,                            &
     &  cool_rate, L_force, L_force_lbc,                                &
     &  zprofile_orog, idl_interp_option, hf,                           &
     &  L_fix_orog_hgt_lbc, orog_hgt_lbc,                               &
     &  L_trivial_trigs, f_plane, ff_plane, r_plane,                    &
     &  idl_bubble_option, idl_bubble_max                               &
     &, idl_bubble_height, idl_bubble_width, idl_bubble_depth           &
     &, idl_bubble_xoffset,idl_bubble_yoffset                           &
     &, L_idl_bubble_saturate,                                          &
     &  L_rotating, L_fixed_lbcs, L_code_test,                          &
     &  L_baroclinic, L_cyclone,                                        &
     &  L_damp, L_geo_for, L_bomex,                                     &
     &  DMPTIM, HDMP, ZDMP,                                             &
     &  u_geo, v_geo
! CRUNTIMC end
! Description: COMDECK containing problem_number
!  for use in setting problem types
!
! Author : T. Davies
! History:
! Version  Date      Comment.
! 5.3      15/11/01  New code

      INTEGER, PARAMETER :: standard=0
      INTEGER, PARAMETER :: monsoon=1
      INTEGER, PARAMETER :: dynamical_core=2
      INTEGER, PARAMETER :: idealised_problem=3
      INTEGER, PARAMETER :: standard_namelist=4
! Description: COMDECK containing surface types
!  for use in idealised problems
!
! Author : T. Davies
! History:
! Version  Date      Comment.
! 5.3      15/11/01  New code

      INTEGER, PARAMETER :: surface_zero=0
      INTEGER, PARAMETER :: surface_ellipse=1
      INTEGER, PARAMETER :: surface_ridge=2
      INTEGER, PARAMETER :: surface_plateau=3
      INTEGER, PARAMETER :: surface_massif=4
      INTEGER, PARAMETER :: surface_mask=5
      INTEGER, PARAMETER :: surface_gauss=6
      INTEGER, PARAMETER :: surface_ridge_series=7
      INTEGER, PARAMETER :: surface_dump=10
! Description: COMDECK containing vertical grid types
!  for use in idealised problems
!
! Author : T. Davies
! History:
! Version  Date      Comment.
! 5.3      15/11/01  New code

      INTEGER, PARAMETER :: vert_regular=1
      INTEGER, PARAMETER :: vert_quadratic_theta=21
      INTEGER, PARAMETER :: vert_bi_quadratic=22
      INTEGER, PARAMETER :: vert_quadratic_uvtheta=23
      INTEGER, PARAMETER :: vert_schar=3
      INTEGER, PARAMETER :: vert_dwd=4
      INTEGER, PARAMETER :: vert_stretch_plus_regular=5
      INTEGER, PARAMETER :: vert_quad_stretch_thin=6
      INTEGER, PARAMETER :: vert_regular_thin=7
      INTEGER, PARAMETER :: vert_dump=10
!+ ---------------------------------------------------------------------
!  Module to specify allowed methods of interpolating from the
!  ancillary file.
!
!  Current Code Owner: J. M. Edwards
!
!  History:
!
!  Version  Date      Comment.
!  5.2      14/11/00  Original code.
!                     (J. M. Edwards)
!
!- ---------------------------------------------------------------------
!
      Integer, parameter :: IO3_3DSPEC = 1
!       Ozone is provided as a full 3D field.
      Integer, parameter :: IO3_2DSPEC = 2
!       Ozone is expanded from a 2D field by direct copying.
      Integer, parameter :: IO3_2DMASSCON = 3
!       Ozone is expanded from a 2D field with conservation of mass.
      Integer, parameter :: IO3_TROP_MAP = 4
!       Ozone mixing ratios are set by mapping each height at each
!       grid-point in the real profile to a height in the ancillary
!       profile and using the mixing ratio there. The mapping is
!       set using the height of the tropopause
      Integer, parameter :: IO3_TROP_MAP_MASSCON = 5
!       Ozone mixing ratios are set as above, but scaled so as to
!       preserve the vertically integrated column ozone
!
! ----------------------------------------------------------------------
! CTIME ----------------------------------------------------
!
!  Purpose: Derived model time/step information including start/end
!           step numbers and frequencies (in steps) of interface field
!           generation, boundary field updating, ancillary field
!           updating; and assimilation start/end times.
!           NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!           Also contains current time/date information, current
!           step number (echoed in history file) and steps-per-group.
!
!  Model            Modification history from model version 3.0:
! version  Date
!
!   3.1   13/02/93  Dimension arrays A_INTERFACE_STEPS/FSTEP/LSTEP
!                   D. Robinson
!   3.3  01/02/94  Add BASIS_TIME_DAYS to BASIS_TIME_SECS for revised
!                  (32-bit portable) model clock calculations. TCJ
!  3.4  13/12/94  Change COMMOM name from CTIME to CTIMED to satisfy
!                 DEC alpha compiler for portability.  N.Farnon.
!  3.5  12/04/95  Stage 1 submodel changes: move to dimensioning
!                 arrays by internal model. R.Rawlins
!  4.4  06/10/97  Data time of IAU dump added. Adam Clayton.
!  4.5  21/08/98  Remove redundant code. D. Robinson.
!  5.1  13/04/00  Instead of saving full IAU data time, save step on
!                 which data time must be reset during an IAU run.
!                 Adam Clayton
!  5.5  17/02/03  Upgrade Wave model from 4.1 to 5.5 D.Holmes-Bell
!
! Programming standard :
!
!  Logical components covered: C0
!
! Project task :
!
! External documentation: Unified Model documentation paper No:
!                         Version:
!
!END -----------------------------------------------------------------

      INTEGER :: I_YEAR               ! Current model time (years)
      INTEGER :: I_MONTH              ! Current model time (months)
      INTEGER :: I_DAY                ! Current model time (days)
      INTEGER :: I_HOUR               ! Current model time (hours)
      INTEGER :: I_MINUTE             ! Current model time (minutes)
      INTEGER :: I_SECOND             ! Current model time (seconds)
      INTEGER :: I_DAY_NUMBER         ! Current model time (day no)
      INTEGER :: PREVIOUS_TIME(7)     ! Model time at previous step
      INTEGER :: IAU_DTResetStep      ! Data time reset step for IAU run

      INTEGER :: BASIS_TIME_DAYS  ! Integral no of days to basis time
      INTEGER :: BASIS_TIME_SECS  ! No of seconds-in-day at basis time
      INTEGER :: O_CLM_FIRSTSTEP  ! First } step for ocean climate
      INTEGER :: O_CLM_LASTSTEP   ! Last  } increments

      LOGICAL :: L_C360DY

! UM6.5MODEL_ANALYSIS_HRS changed to REAL - 
!   requires FORECAST_HRS and DATA_MINUS_BASIS_HRS to REAL also 
      REAL    :: FORECAST_HRS     ! Hours since Data Time (ie T+nn)
      REAL    :: DATA_MINUS_BASIS_HRS ! Data time - basis time (hours)

      COMMON /CTIMED/ I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,    &
     &  I_DAY_NUMBER,PREVIOUS_TIME,                                     &
     &  BASIS_TIME_DAYS,BASIS_TIME_SECS,                                &
     &  FORECAST_HRS,DATA_MINUS_BASIS_HRS,                              &
     &  IAU_DTResetStep,                                                &
     &  O_CLM_FIRSTSTEP,   O_CLM_LASTSTEP, L_C360DY

      INTEGER :: STEPim(INTERNAL_ID_MAX)  ! Step no since basis time
      INTEGER :: GROUPim(INTERNAL_ID_MAX) ! Number of steps per group

      ! Finish step number this run
      INTEGER :: TARGET_END_STEPim(INTERNAL_ID_MAX)

      REAL :: SECS_PER_STEPim(INTERNAL_ID_MAX) ! Timestep length in secs

      ! Frequency of interface field generation in steps
      INTEGER :: INTERFACE_STEPSim(MAX_N_INTF,INTERNAL_ID_MAX)

      ! Start steps for interface field generation
      INTEGER :: INTERFACE_FSTEPim(MAX_N_INTF,INTERNAL_ID_MAX)

      ! End steps for interface field generation
      INTEGER :: INTERFACE_LSTEPim(MAX_N_INTF,INTERNAL_ID_MAX)

      ! Frequency of  updating boundary fields in steps
      INTEGER :: BOUNDARY_STEPSim(INTERNAL_ID_MAX)

      ! No of steps from boundary data prior to basis time to model
      ! basis time
      INTEGER :: BNDARY_OFFSETim(INTERNAL_ID_MAX)

      ! Lowest frequency for updating of ancillary fields in steps
      INTEGER :: ANCILLARY_STEPSim(INTERNAL_ID_MAX)

      ! Start steps for assimilation
      INTEGER :: ASSIM_FIRSTSTEPim(INTERNAL_ID_MAX)

      ! Number of assimilation steps to analysis
      INTEGER :: ASSIM_STEPSim(INTERNAL_ID_MAX)

      ! Number of assimilation steps after analysis
      INTEGER :: ASSIM_EXTRASTEPSim(INTERNAL_ID_MAX)

      COMMON/CTIMEE/                                                    &
     &  STEPim,GROUPim,TARGET_END_STEPim,INTERFACE_STEPSim,             &
     &  INTERFACE_FSTEPim,INTERFACE_LSTEPim,BOUNDARY_STEPSim,           &
     &  BNDARY_OFFSETim,ANCILLARY_STEPSim,ASSIM_FIRSTSTEPim,            &
     &  ASSIM_STEPSim,ASSIM_EXTRASTEPSim,SECS_PER_STEPim

! CTIME end
!
! Description:
!   This comdeck declares an integer variable 'ModelType' whose value
!   determines whether a model run is global, limited area or zonal.
!   The values of ModelType associated with each of the run types are
!   defined by integer parameters which are also declared below.
!    ModelType is set in subroutine SETLOGIC.
!
! Current Code Owner: S.J.Swarbrick
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   3.4    29/9/94  Original code. S.J.Swarbrick
!
      ! Value used to represent the global model configuration
      INTEGER,PARAMETER:: GlobalModel=1

      ! Value used to represent the limited area model configuration
      INTEGER,PARAMETER:: LimitedAreaModel=2

      ! Value used to represent the 'periodic in x' model config
      ! INTEGER,PARAMETER:: ZonalModel=2

! Global scalars:
      INTEGER     ModelType  ! Integer switch which is equated to one
                             ! of the above parameters in a model run,
                             ! and so determines the configuration

! COMMON blocks:
      COMMON /RunType/ ModelType

! C_GLOBAL end
!*L -----------------COMDECK PHYSCONS----------------------------------
!
!  Purpose : contains physical constants required by the whole of the
!            model. It is made up of individual COMDECKS for sets of
!            of related constants, each routine can access one or
!            several of these COMDECKS seperately
!  System Component : F07
!  System task : Z
! END
!*----------------------------------------------------------------------
!*L------------------COMDECK C_A----------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Replace variable A by more meaningful name for
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Convert to Fixed/Free format. P. Selwood

      ! Mean radius of Earth in m.
      Real, Parameter  :: Earth_Radius = 6371000.

!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!*L------------------COMDECK C_OMEGA------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable two_omega for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.3   12/10/01  two_omega initialised in SETCON     Terry Davies
!OMEGA IS MAGNITUDE OF EARTH'S ANGULAR VELOCITY
      REAL OMEGA                                                        &
     &,two_omega

       Common/Omega/Omega, two_omega
!  Angular speed of Earth's rotation Omega to be initialised in SETCON
!*----------------------------------------------------------------------
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
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
! C_KT_FT start

      REAL,PARAMETER:: KT2MS=1852.0/3600.0 ! Knots to m/s conversion
      REAL,PARAMETER:: FT2M =0.3048        ! Feet to meters conversion

! C_KT_FT end
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
!*L------------------COMDECK C_LAPSE ----------------------------------
      Real, Parameter :: lapse      = 0.0065  ! Near surface lapse rate
      Real, Parameter :: lapse_trop = 0.002   ! Tropopause lapse rate
!*---------------------------------------------------------------------
!
! This Comdeck declares and stores the logical and integer
! variables used in the time-step control for writing general
! data.
!
! Code owner: S.J.Swarbrick
!
! Switch which activates output of arrays for Peer utility
      LOGICAL L_PEER

! Switches which activate dump writing
      LOGICAL                                                           &
     &  L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                  &
     &  ,L_WRIT_INIT

! Timesteps for dump writing
      INTEGER                                                           &
     &  T_WRITD1_START                                                  &
                              ! First timestep
     &  ,T_WRITD1_END                                                   &
                              ! Last timestep
     &  ,T_WRITD1_INT         ! Timestep interval between dumps

      INTEGER                                                           &
     &  PEER_VN                  ! Version of PEER utility

      NAMELIST/NLSTWRITEDATA/                                           &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT

      COMMON/WRITEDATA/                                                 &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT
!
! COMDECK PPXLOOK
! Description:
!
!   Declares ppxref look-up arrays used by the UM and associated
!    arrays and parameters.
!   Comdecks CSUBMODL,CPPXREF must be *CALLed before this
!    comdeck
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       May. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.0       Dec. 95   Replace dynamic dim ppxRecs with
!                     NUM_DIAG_MAX in PPXC   N. Farnon
! 4.1       July 96   *CALL VERSION introduced - NUM_DIAG_MAX made
!                      equal to NDIAGP.
!                     NUM_USR_DIAG_MAX increased from 200 to 300
!                      (just in case).
! 4.4       03/11/97  Removed MKPPXRF *DEF references. K Rogers
! 4.4       04/11/97  Changed -RECON def line to allow for other small
!                     execs which had used the RECON def. K Rogers
!
! Declarations:

! Global parameters:
! VERSION STASH parameter definitions
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Mar. 95   Original code.  S.J.Swarbrick
! 4.0                                 S.J.Swarbrick
! 4.1       Apr. 96   Rationalise MDI  S.J.Swarbrick
!  4.1  29/02/96  Increase OUTFILE_E.  RTHBarnes.
!  4.2  27/11/96  mpp code : Increase NELEMP   P.Burton
!  4.3  04/06/97  Increase NELEMP for D1 addressing S.D.Mullerworth
!  4.4  04/06/97  Increase NELEMP for sampling offset. S.D.Mullerworth
!  4.5  28/01/98  Increade NELEMP for mpp code.   P.Burton
!  4.5  18/09/98  Modify name of VERSION common block to stop potential
!                 clashes with Fortran variable names          P.Burton
!  4.5  30/09/98  Increase NRECDP from 600 to 800. D. Robinson.
!  5.2  29/01/01  OUTFILE_E changed. Adam Clayton
!  5.5  20/02/03  Increased size of STASH_SET.  P.Dando
!  6.1  03/08/04  Increase size of NPSLEVP and NPSLISTP
!                 (Pseudo Levels)  Anthony A. Dickinson
!  6.1  04/08/04  Increase size of NDIAGP W Roseblade.
!  6.2  31/03/06  Increase size of NDIAGP again.
!                 R Sempers (frpz)
!  6.2  06/04/06  Increased size of OUTFILE_E   T. Edwards
!  6.2  03/02/06  Increase NRECDP to 1500. T Johns
!
      ! Max. no. of STASH sections  per internal model (44 in practice)
      INTEGER,PARAMETER :: NSECTP=99
      ! Max. no. of STASH items per section
      INTEGER,PARAMETER :: NITEMP=999
      ! Max. no. of STASH list records (prognostic + diagnostic)
      INTEGER,PARAMETER :: NRECDP=1500
      ! Max. no. of output times tables in STASHC
      INTEGER,PARAMETER :: NTIMEP=100
      ! Max. no. of time profiles in STASHC
      INTEGER,PARAMETER :: NPROFTP=100
      ! Max. no. of domain profiles/levels lists in STASHC (used for
      ! both)
      INTEGER,PARAMETER :: NPROFDP=100
      ! Max. total no. of time series in STASHC
      INTEGER,PARAMETER :: NTimSerP=1500
      ! Max. no. time series per domain profile
      INTEGER,PARAMETER :: tsdp=250
      ! Max. no. of useage profiles in STASHC
      INTEGER,PARAMETER :: NPROFUP=40
      ! Max. no. of levels in a levels list
      INTEGER,PARAMETER :: NLEVP=50
      ! Max. no. of pseudo levels in a  pseudo levels list
      INTEGER,PARAMETER :: NPSLEVP=100
      ! Max. no. of pseudo levels lists in STASHC
      INTEGER,PARAMETER :: NPSLISTP=100
      ! Max. no. non-blank records in PPXREF file
      INTEGER,PARAMETER :: NDIAGP=2600
      INTEGER,PARAMETER :: NDIAGPM=NRECDP  ! Same as NRECDP
      INTEGER,PARAMETER :: NELEMP=33
      INTEGER,PARAMETER :: NLEVP_S=NLEVP*6+1
      INTEGER,PARAMETER :: NLEVLSTSP=NPROFDP
      INTEGER,PARAMETER :: NMEANP=4  ! No. of meaning periods
      ! OUTFILE_S, OUTFILE_L and OUTFILE_E must be consistent with
      ! NUNITS and NUNITS_LEN in file CHSUNITS.
      ! Ranges of output file numbers
      INTEGER,PARAMETER :: OUTFILE_S=20
      INTEGER,PARAMETER :: OUTFILE_E=161
      INTEGER,PARAMETER :: OUTFILE_L=OUTFILE_E-OUTFILE_S+1
!Global scalar:
      CHARACTER(LEN=80) :: STASH_SET     !Names of stasets files
!Common block:
      COMMON/common_VERSION/ STASH_SET
! VERSION end
! No. of STASH items per section
      INTEGER      PPXREF_ITEMS
        PARAMETER (PPXREF_ITEMS    =NITEMP)
! No. of STASH sections per internal model
      INTEGER      PPXREF_SECTIONS
        PARAMETER (PPXREF_SECTIONS =NSECTP-55)
! Max. number of non-null records in ppxref file (>1200)
      INTEGER      NUM_DIAG_MAX
        PARAMETER (NUM_DIAG_MAX    =NDIAGP)
! Max. number of user-defined ppxref records allowed
      INTEGER      NUM_USR_DIAG_MAX
        PARAMETER (NUM_USR_DIAG_MAX=450)

! No. of ppxref records read into PPXI,PPXC (for dyn. allocation)
      INTEGER      ppxRecs

! Global arrays:
! ppxref look-up arrays
      INTEGER   PPXI(ppxRecs,PPXREF_CODELEN)
      CHARACTER PPXC(NUM_DIAG_MAX,PPXREF_CHARLEN)
! Arrays for temporary storage of user-ppxref records -
!   used to transfer these records from STASH_PROC into U_MODEL
      INTEGER   PPXI_U(NUM_USR_DIAG_MAX,PPXREF_CODELEN)
      CHARACTER PPXC_U(NUM_USR_DIAG_MAX,PPXREF_CHARLEN)
! Array of flags to indicate origin of ppxref record
! 'P' for ppxref file; 'U' for user-stash master file
      CHARACTER OriginFlag(NUM_DIAG_MAX)
! Array of indices to identify which ppxref record corresponds to
!   any given row of PPXI, PPXC
      INTEGER   RowIndex(NUM_DIAG_MAX)
! Pointer array for PPXI, PPXC arrays
      INTEGER PPXPTR                                                    &
     & (N_INTERNAL_MODEL    ,0:PPXREF_SECTIONS ,PPXREF_ITEMS)

! Common block:
      COMMON/PPX_INT/ RowIndex,PPXI_U
      COMMON/PPX_CHA/ OriginFlag,PPXC_U
! - End --------------------------------------------------------------
! CPRINTST defines print status for standard output messages

      ! Minimum output, only essential messages
      INTEGER,PARAMETER :: PrStatus_Min    = 1

      ! Normal informative messages + warnings
      INTEGER,PARAMETER :: PrStatus_Normal = 2

      ! Operational status, all informative messages
      INTEGER,PARAMETER :: PrStatus_Oper   = 3

      ! All informative + extra diagnostic messages
      INTEGER,PARAMETER :: PrStatus_Diag   = 4

      INTEGER PrintStatus ! Control volume of standard output messages
      COMMON/PrintSt/PrintStatus

! CPRINTST end
!C_DUST_NDIV.............................................................
! Description: Contains parameters for mineral dust code
! Current Code Owner: Stephanie Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  5.5      12/02/03  Original Code.   Stephanie Woodward
!
! Declarations:
!
      INTEGER NDIV        ! number of particle size divisions
      PARAMETER (NDIV = 6)
!.....................................................................
! ----------------------- header file: CTFilt  -----------------------
!
! Description:
!
!   Parameters and variables for Incremental Analysis Update (IAU) and
!   Temporal Digital Filtering (TDF) schemes.
!
!
! Current Code Owner: Adam Clayton.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      INTEGER :: TDF_unit,                                              &
                             ! Unit number for TDF output dump.
     &           IAU_unit    ! Unit number for IAU increment.

      PARAMETER (TDF_unit = 107)
      PARAMETER (IAU_unit = 108)

      INTEGER :: LeaveAlone,                                            &
                             ! \  Options for filtering
     &           DirectFilt,                                            &
                             !  > of advected wind
     &           CopyFiltUVW ! /  components.

      PARAMETER (LeaveAlone  = 1)
      PARAMETER (DirectFilt  = 2)
      PARAMETER (CopyFiltUVW = 3)

      REAL    :: q_CC_tol,                                              &
                             ! qCL/qCF tolerance for cloud clearing.
     &           Weight_LL   ! Lower limit for absolute value of
                             ! filter weight.

      PARAMETER (q_CC_tol  = 1.0E-12  )
      PARAMETER (Weight_LL = 0.0000001)

      REAL    :: q_min       ! Minimum value allowed for q after
                             ! addition of q increments.

      REAL    :: oz_min      ! Minimum value allowed for ozone after
                             ! addition of ozone increments

      INTEGER :: MaxNumWeights ! Maximum number of filter weights.

      PARAMETER (MaxNumWeights = 1000)

      REAL    :: IAU_Weights(MaxNumWeights),                            &
                                             ! IAU filter weights array.
     &           TDF_Weights(MaxNumWeights)  ! TDF filter weights array.

      INTEGER :: IAU_FixHd(LEN_FIXHD) ! Fixed-length header of
                                      ! IAU increment dump.

      INTEGER :: D1_TFilt_len,                                          &
                                 ! D1 length for unpacked fields
     &           D1_IAU_k4_len,                                         &
                                 ! D1 length for packed IAU fields
     &           IAU_Len1Lookup,                                        &
                                 ! \ IAU increment
     &           IAU_Len2Lookup  ! / lookup dimensions

      LOGICAL :: L_Pack_D1_IAU   ! Hold IAU fields in packed form?

      INTEGER :: IAU_NumFldCodes
      PARAMETER (IAU_NumFldCodes = 13)

      ! Codes of fields that may be read in from the IAU file:
      INTEGER :: IAU_FldCodes(IAU_NumFldCodes)
      PARAMETER (IAU_FldCodes = (/ 2,   3,   4,   10,  12, 90,          &
     &                             150, 253, 254, 255, 407, 18001,480 /))

      ! Local lengths of fields (0 if field not required):
      INTEGER :: IAU_LocFldLens(IAU_NumFldCodes)

      ! Field descriptions:
      CHARACTER(7) :: IAU_FldDescs(IAU_NumFldCodes)
      PARAMETER  (IAU_FldDescs = (/ 'u      ', 'v      ', 'theta  ',    &
     &                              'q      ', 'qCF    ', 'aerosol',    &
     &                              'w      ', 'rho    ', 'qCL    ',    &
     &                              'exner  ', 'p      ', 'qT     ',     &
     &                              'ozone  ' /))

  ! IAU namelist variables:
  ! -----------------------

      LOGICAL L_IAU               ! Activate IAU scheme?

      INTEGER IAU_StartMin,                                             &
                                  ! Start minute of filtering period.
     &        IAU_EndMin,                                               &
                                  ! End   minute of filtering period.
     &        IAU_ApexMin         ! Apex minute for triangular filter.

      REAL    IAU_Cutoff_period,                                        &
                                  ! Filter cut-off period in hours.
     &        IAU_SBE_period      ! Stop band edge period in hours.

      LOGICAL L_IAU_CalcExnerIncs ! Use p increments to calculate
                                  ! exner increments?

      LOGICAL L_IAU_CalcThetaIncs ! Calculate theta increments using
                                  ! exner and q increments?

      LOGICAL L_IAU_CalcRhoIncs   ! Calculate rho increments using
                                  ! exner, theta and (possibly) q
                                  ! increments?

      LOGICAL L_IAU_IncTStar      ! If set, add level-one temperature
                                  ! increments to surface temperature
                                  ! and top-level soil temperature.

      LOGICAL L_IAU_ResetPoles    ! If set, reset polar rows of
                                  ! relevant increment fields to
                                  ! their mean values.

      LOGICAL L_IAU_RemoveSS      ! Remove supersaturation wrt water?

      LOGICAL L_IAU_CallStratQ    ! Reset stratospheric humidities at
                                  ! end of IAU insertion period?

      LOGICAL L_IAU_Diags         ! If set, write out IAU diagnostics.

      LOGICAL L_IAU_LowMem        ! If set, activate low memory (but
                                  ! high IO) version of the IAU code.

      LOGICAL L_IAU_DumpTS0State  ! If set, write out model state
                                  ! immediately after timestep-zero
                                  ! call to TFilt_cntl.

      LOGICAL L_IAU_CalcCloudIncs ! If set, calculate q, qcl & Cl
                                  ! from q or qT increments.

      LOGICAL L_IAU_IncrementIce  ! If this and L_IAU_CalcCloudIncs,
                                  ! calculate qcf & Cf

      LOGICAL L_IAU_ScaleCloud    ! If set, scale qcl, qcf, Cl & Cf
                                  ! increments to be in physical bounds
                                  
      LOGICAL L_IAU_UPPER_THETA   ! If set, then constrain the upper theta 
                                  ! increments.

      LOGICAL L_IAU_SetOzoneMin  ! If set, reset ozone to oz_min        
                                 ! in IAU if ozone was negative

      CHARACTER*10                                                      &
     &        IAU_FilterType      ! Filter type.

      REAL    IAU_LL_strat_pv,                                          &
                                  ! Lower-limit for strat ABS(PV).
     &        IAU_UL_strat_p,                                           &
                                  ! Upper-limit for strat pressure.
     &        IAU_LL_trop_p       ! Lower-limit for trop  pressure.

  ! TDF namelist variables:
  ! -----------------------

      LOGICAL L_TDF               ! Activate TDF scheme?

      INTEGER TDF_StartMin,                                             &
                                  ! Start minute of filtering period.
     &        TDF_EndMin,                                               &
                                  ! End   minute of filtering period.
     &        TDF_ApexMin         ! Apex minute for triangular filter.

      REAL    TDF_Cutoff_period,                                        &
                                  ! Filter cut-off period in hours.
     &        TDF_SBE_period      ! Stop band edge period in hours.

      LOGICAL L_TDF_FilterQ,                                            &
                                  ! Filter q?
     &        L_TDF_FilterQCL,                                          &
                                  ! Filter qCL?
     &        L_TDF_FilterQCF     ! Filter qCF?

      LOGICAL L_TDF_ModifyCloud   ! If set, modify cloud variables in
                                  ! TDF dump so that it becomes suitable
                                  ! for starting forecasts including
                                  ! physics.

      INTEGER TDF_AdvWindOpt      ! Filtering option for advected winds.

      LOGICAL L_TDF_CallStratQ    ! Reset stratospheric humidities in
                                  ! TDF dump?

      CHARACTER*10                                                      &
     &        TDF_FilterType      ! Filter type.


      COMMON / CTFilt /                                                 &
     & q_min, oz_min, IAU_Weights, TDF_Weights, IAU_FixHd,              &
     & D1_TFilt_len, D1_IAU_k4_len, IAU_Len1Lookup, IAU_Len2Lookup,     &
     & L_Pack_D1_IAU, IAU_LocFldLens,                                   &
     & L_IAU, IAU_StartMin, IAU_EndMin,                                 &
     & IAU_ApexMin, IAU_Cutoff_period, IAU_SBE_period,                  &
     & L_IAU_CalcCloudIncs, L_IAU_IncrementIce, L_IAU_ScaleCloud,       &
     & L_IAU_CalcExnerIncs, L_IAU_CalcThetaIncs, L_IAU_CalcRhoIncs,     &
     & L_IAU_IncTStar, L_IAU_ResetPoles, L_IAU_RemoveSS,                &
     & L_IAU_CallStratQ, L_IAU_Diags, L_IAU_LowMem, L_IAU_DumpTS0State, &
     & L_IAU_UPPER_THETA,                                               &
     & IAU_LL_strat_pv, IAU_UL_strat_p, IAU_LL_trop_p,                  &
     & L_TDF, TDF_StartMin, TDF_EndMin,                                 &
     & TDF_ApexMin, TDF_Cutoff_period, TDF_SBE_period,                  &
     & L_TDF_FilterQ, L_TDF_FilterQCL, L_TDF_FilterQCF,                 &
     & L_TDF_ModifyCloud, TDF_AdvWindOpt,                               &
     & L_TDF_CallStratQ, L_IAU_SetOzoneMin,                             &
     ! Character variables at the end of common block
     & IAU_FilterType, TDF_FilterType

      NAMELIST / RUN_TFilt /                                            &
     & L_IAU, IAU_StartMin, IAU_EndMin,                                 &
     & IAU_ApexMin, IAU_Cutoff_period, IAU_SBE_period,                  &
     & L_IAU_CalcCloudIncs, L_IAU_IncrementIce, L_IAU_ScaleCloud,       &
     & L_IAU_CalcExnerIncs, L_IAU_CalcThetaIncs, L_IAU_CalcRhoIncs,     &
     & L_IAU_IncTStar, L_IAU_ResetPoles, L_IAU_RemoveSS,                &
     & L_IAU_CallStratQ, L_IAU_Diags, L_IAU_LowMem, L_IAU_DumpTS0State, &
     & L_IAU_UPPER_THETA,                                               &
     & IAU_LL_strat_pv, IAU_UL_strat_p, IAU_LL_trop_p,                  &
     & L_TDF, TDF_StartMin, TDF_EndMin,                                 &
     & TDF_ApexMin, TDF_Cutoff_period, TDF_SBE_period,                  &
     & L_TDF_FilterQ, L_TDF_FilterQCL, L_TDF_FilterQCF,                 &
     & L_TDF_ModifyCloud, TDF_AdvWindOpt,                               &
     & L_TDF_CallStratQ, L_IAU_SetOzoneMin,                             &
     & IAU_FilterType, TDF_FilterType
! Start C_KINDS

! Description:
!   Contains parameters defining kinds for 32 and 64 integers
!   and reals.
!
! Current Code Owner: Paul Selwood
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   6.1  10/03/04   Original code. Paul Selwood

! Parameters for 32 and 64 bit kinds

! Precision and range for 64 bit real
      Integer, Parameter :: prec64  = 15
      Integer, Parameter :: range64 = 307

! Precision and range for 32 bit real
      Integer, Parameter :: prec32  = 6
      Integer, Parameter :: range32 = 37

! Range for integers
      Integer, Parameter :: irange64=15
      Integer, Parameter :: irange32=9

! Kind for 64 bit real
      Integer, Parameter :: real64  = selected_real_kind(prec64,range64)
! Kind for 32 bit real
      Integer, Parameter :: real32  = selected_real_kind(prec32,range32)
! Kind for 64 bit integer
      Integer, Parameter :: integer64 = selected_int_kind(irange64)
! Kind for 32 bit integer
      Integer, Parameter :: integer32 = selected_int_kind(irange32)

! Kinds for 64 and 32 bit logicals. Note that there is no
! "selected_logical_kind", but using the equivalent integer kind is a
! workaround that works on every platform we've tested.
      Integer, Parameter :: logical64 = integer64
      Integer, Parameter :: logical32 = integer32

! End C_KINDS

! Variables for FLUME
      integer :: tr_size ! size of a single tracer

! Subroutine arguments:
      INTEGER :: IAU_lookup(IAU_Len1Lookup, IAU_Len2Lookup)

      REAL(KIND=real32) :: D1_IAU_k4(D1_IAU_k4_len) ! Packed IAU fields
      REAL         :: D1_TFilt (D1_TFilt_len ) ! Unpacked IAU/TDF fields
      
! 3-D fields of species to be passed down to radiation       
      INTEGER, INTENT(IN) :: ngrgas 
      INTEGER, INTENT(IN) :: grgas_addr(ngrgas)

! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Atm_Step')

      LOGICAL gather   ! Convert to sea_ice points within BL code
      PARAMETER ( gather = .True. ) ! (was L_COMPRESS_SEAICE)

! Local scalars:
      Integer                                                           &
     &  i, j, k, l                                                      &
                          ! loop counters
     &, CycleNo                                                         &
                   ! Number of cycles (iterations) for iterative SISL.
     &, ij                                                              &
                          ! 2D array index for offx/y variables
     &, lbc_pt                                                          &
                          ! pointer to LBC array
     &, ij_u                                                            &
                          ! 2D array index for U offx/y variables
     &, ij_v                                                            &
                          ! 2D array index for V offx/y variables
     &, ji                                                              &
                          ! 2D array index for halo_i/j variables
     &, gi                                                              &
                          ! derived from loop counters
     &, lambda_start                                                    &
                    ! pointer for start of lambda_p/lambda_u on this pe

     &, first_constant_r_rho_level_m1                                   &
                                      ! value used to dimension
!                                     ! arrays, max of (1 and
!                                     ! first_constant_r_rho_level)
     &, icount                                                          &
                ! diagnostic count
     &, i_start, i_stop                                                 &
     &, j_start, j_stop                                                 &
     &, j_begin, j_end                                                  &
     &, j_ti, j_ice_fract, j_ice_thick

      Integer                                                           &
     &  n_Y_arrays                                                      &
                         ! = 1 for global, 3 for LAM
     &, n_Yw_arrays                                                     &
                         ! = 1 for global, 2 for LAM
     &, n_Yd_arrays                                                     &
                         ! = 1 for global, 3 for LAM
     &, n_Ydw_arrays     ! = 1 for global, 2 for LAM
!
      Integer :: nd_o3  ! Total size of ozone array supplied
      Integer :: nd_stochem  ! size of stochem ozone array

      Integer                                                           &
     &  check_bottom_levels ! used in interpolation code, and is
!                           ! the number of levels to check to see
!                           ! if the departure point lies inside the
!                           ! orography.

      INTEGER                                                           &
     &  timestep_number                                                 &
                             ! no. of atmos timesteps since basis time
     &,first_constant_r_rho_level                                       &
                                  ! 1st rho level on which r is constant
     &, info       ! icode return from UM_FORT_FLUSH

      INTEGER                                                           &
     &  DIM_CS1                                                         &
                         ! soil C dimension: 1 for single, 4 for RothC
     &, DIM_CS2          ! soil C dimension: 1 for single,
                         !                        LAND_FIELD for RothC
      REAL                                                              &
     & constant                                                         &
     &,h_print                                                          &
                         ! needed for printing idealised orography
     &,timestep                                                         &
                            ! atmosphere model timestep
     &,radiation_timestep                                               &
                            ! timestep of radiation scheme
     &,radiation_tstep_diag                                             &
                              ! timestep of fast radiation (3C)
     &,radiation_tstep_prog                                             &
                              ! timestep of slow radiation (3C)
     &,chemistry_timestep             ! must be  <=  model timestep

      REAL                                                              &
     & pos_timestep                                                     &
                    ! = +timestep.
     &,neg_timestep ! = -timestep.

      LOGICAL                                                           &
     & l_rad_step                                                       &
                            ! :T :activate radiation this timestep
     &,L_Rad_step_diag                                                  &
                         ! T :activate fast radiation this timestep (3C)
     &,L_Rad_step_prog                                                  &
                         ! T :activate slow radiation this timestep (3C)
     &,L_print_L2norms                                                  &
                         ! diagnostic printing of l2norms
     &,L_print_L2helm                                                   &
                         ! diagnostic printing of l2norms in solver
     &,L_Tracer          ! T if *any* tracer variables present

! Code to do with tropopause diagnostics from O3_to_3D
! Declare logicals for Stash Diagnostics used in call to O3_to_3D
      Logical :: L_O3_trop_level   ! STASH code 2,280
      Logical :: L_O3_trop_height  ! STASH code 2,281
      Logical :: L_T_trop_level    ! STASH code 2,282
      Logical :: L_T_trop_height   ! STASH code 2,283

! Declare the tropopause variables output for O3_to_3D as allocatable
      Real, DIMENSION (:,:), ALLOCATABLE::                              &
     & O3_trop_level, T_trop_level, O3_trop_height, T_trop_height

! END OF Code to do with tropopause diagnostics from O3_to_3D

      INTEGER                                                           &
     & LAND_PTS_TRIF                                                    &
                     !\ For dimensioning variables in NI_bl_ctl
     &,NPFT_TRIF                                                        &
                     !/ depending on whether TRIFFID is in use
     &,CO2_DIM_LEN                                                      &
                     !\ For dimension 3-D CO2 field to be passed
     &,CO2_DIM_ROW                                                      &
                     !/ to NI_bl_ctl
     &,CO2_DIM_LEV   !/ and NI_rad_ctl

      INTEGER :: Sec         ! Second
      LOGICAL :: L_CallTFilt ! Call TFilt_cntl?

! Variables required for call to SET_LATERAL_BOUNDARIES
      INTEGER                                                           &
     &  lbc_size           ! size of a single level of LBC

      LOGICAL                                                           &
     &  L_do_halos                                                      &
                          ! update the halos?
     &, L_do_boundaries
                          ! update the boundaries?


! Variables required for ice category selection
      REAL, POINTER :: p_ti(:,:,:)
      REAL, POINTER :: p_ice_fract(:,:,:)
      REAL, POINTER :: p_ice_thick(:,:,:)

! soil carbon content & accumulated soil respiration 
! (target varies according to L_TRIFFID)
      REAL, POINTER :: CS(:), RSA(:)

! Array dimensions for sea-salt aerosol
      INTEGER                                                           &
     &  salt_dim1                                                       &
     &, salt_dim2                                                       &
     &, salt_dim3

! Array dimensions for Aero_Ctl
      INTEGER                                                           &
     & aero_dim1                                                        &
     &,aero_dim2                                                        &
     &,aero_dim3

! CABLE declarations
! C_SOILH start
      ! No. of soil layers (must = NSOIL).
      REAL,PARAMETER:: PSOIL=4

      ! Tunable characteristic freq (rad/s)
      REAL,PARAMETER:: OMEGA1=3.55088E-4

      ! Density of lying snow (kg per m**3)
      REAL,PARAMETER:: RHO_SNOW=250.0

      ! Depth of `effective' snow surface layer (m)
      REAL,PARAMETER:: DEFF_SNOW=0.1

      ! Thermal conductivity of lying snow (Watts per m per K).
      REAL,PARAMETER:: SNOW_HCON=0.265

      ! Thermal capacity of lying snow (J/K/m3)
      REAL,PARAMETER:: SNOW_HCAP=0.63E6

! C_SOILH end
      INTEGER                               &
     &  TILE_INDEX(land_points,NTYPE)       &! OUT Index of tile points     
     &, TILE_PTS(NTYPE)                     &! OUT Number of tile points
     &, day                                 &
     &, total_nsteps                        &! Total number of steps in run
!     &, istep_cur                           &
     &, SOIL_TYPE(row_length,rows)          &
     &, VEG_TYPE(row_length,rows)
      INTEGER, save :: istep_cur=0          ! current time step  
          
      LOGICAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: L_TILE_PTS
!      LOGICAL l_cable
!*** note in ni-rad_ctl & glue_rad alb_tile(land_field,ntiles) dim ***

! CABLE prognostic variables
! Declared as allocable here so they can be saved. Eventually they'll all
! move to the D1 array
      Real, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::                      &
     &  SNOW_COND                                                       &
     &, STHU_TILE
      Real, DIMENSION(:,:), ALLOCATABLE, SAVE ::                        &
     &  T_SURF_TILE

!code differences:
      ! Some of these should eventually be promoted to proper prognostic
      ! variables. Just save them for now - MRD
      ! Save doesn't work because they're automatic variables.
!!$      save SNAGE_TILE, RTSOIL_TILE, GFLUX_TILE, SGFLUX_TILE, SNOW_DEPTH3L    &
!!$     , SNOW_MASS3L, SNOW_COND, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SMCL_TILE &
!!$     , STHU_TILE, STHF_TILE, TSOIL_TILE, T_SURF_TILE, HCONS, alb_tile        &
!!$     , TOT_ALB, CH_CAB, CD_CAB, U_S_CAB, LAND_ALBEDO_CABLE
!code differences:
!      Call RUN_INIT(                                ...                 &
!     &  ,l_use_seasalt_direct, l_use_biogenic, l_murk_rad, l_use_aod    &
!atmos_physics2 needs:

      Real, Dimension(:,:,:), ALLOCATABLE, SAVE ::                            & 
      ALB_TILE,SURF_DOWN_SW,LAND_ALBEDO_CABLE
 
      Real, Dimension(:,:), ALLOCATABLE, SAVE ::                              & 
      FTL_TILE_CAB,LE_TILE_CAB,TSTAR_TILE_CAB,SMCL_CAB,TSOIL_CAB,RADNET_TILE, & 
      SW_DOWN,RTSOIL_TILE,GFLUX_TILE,SGFLUX_TILE,TOT_ALB,CD,CHX,TRANSP_TILE,  &
      WB_LAKE
 
      Real, Dimension(:), ALLOCATABLE, SAVE ::                                & 
      FTL_CAB,LE_CAB,TSTAR_CAB,USTAR_CAB,SURF_HTF_CAB,CH_CAB,CD_CAB,U_S_CAB,  & 
      HCONS,ALBSOIL,TOT_WBLAKE,TOT_SUBRUN

      Real, allocatable, save ::                                        &
     & wblake_ratio    

      Real                                                              &
                lat(row_length, rows),long(row_length, rows),time_sec
                                          ! Lat. of gridpoint chosen
                                          ! Long. of gridpoint chosen
                                          ! actual time of day in secs.
! CABLE end ----------------------------------------------------

! Local Arrays
      Real                                                              &
     &  theta_star(1-offx:row_length+offx,                              &
     &               1-offy:rows+offy, model_levels)                    &
     &, q_star(1-offx:row_length+offx,                                  &
     &           1-offy:rows+offy, wet_levels)                          &
     &, qcl_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, qcf_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, cf_star(1-offx:row_length+offx,                                 &
     &           1-offy:rows+offy, wet_levels)                          &
     &, cfl_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, cff_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, exner_star(1-offx:row_length+offx,                              &
     &             1-offy:rows+offy, wet_levels)                        &
     &, temp_ozone(1-offx:row_length+offx,                              &
     &             1-offy:rows+offy, model_levels)                      &
     &, frac_control(land_field,ntype)   !Forcing for land surface (3C)

      Real, DIMENSION (:,:,:), ALLOCATABLE ::                           &
        ! Local Arrays to store microphysics fields
     &  qcf2_star, qrain_star, qgraup_star

! Declare allocatable arrays for passing cloud fractions
! to LS_ACF_Brooks
      Real, DIMENSION (:,:,:), ALLOCATABLE::                            &
     & cf_bulk_nohalo, cf_liquid_nohalo, cf_frozen_nohalo

      Real                                                              &
     &  R_u(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &        model_levels)                                             &
     &, R_v(1-offx:row_length+offx, 1-offy:n_rows+offy,                 &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)

      Real                                                              &
     &  rho_n (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)

      Real                                                              &
     &   biogenic(row_length, rows, model_levels)

! Local arrays for using the aerosol climatology for NWP

! arcl_dim.h
!
! Maximum dimensions for the aerosol climatology for NWP
!

      integer, parameter :: NPD_ARCL_SPECIES = 7
      integer, parameter :: NPD_ARCL_COMPNTS = 20

! end of arcl_dim.h
      
      ! Internal model switches
      Logical L_USE_ARCL(NPD_ARCL_SPECIES)
      
      ! Internal array of mass-mixing ratios
      Real, dimension(:,:,:,:), allocatable :: arcl
      
      ! Number of requested species within the climatology
      Integer n_arcl_species
      
      ! Corresponding number of requested components
      Integer n_arcl_compnts
      
      ! Array index of each component
      Integer i_arcl_compnts(NPD_ARCL_COMPNTS)
      
      Integer cloud_tol

      Real                                                              &
     &  mag_vector_np (model_levels)                                    &
     &, dir_vector_np (model_levels)                                    &
     &, mag_vector_sp (model_levels)                                    &
     &, dir_vector_sp (model_levels)                                    &
     &, lambda_a (row_length) ! delta_lambda term for polar wind

      REAL, DIMENSION (:,:,:), ALLOCATABLE::                            &
     &  inc_u, inc_v, inc_w, inc_t, inc_rho                             &
     &, inc_q, inc_qcl, inc_qcf, inc_cf, inc_cfl, inc_cff

      Logical                                                           &
     &  L_do_inc_vels                                                   &
     &, L_do_inc_rho               ! flag for rho incr diagnostic

! LAM LBC tendency
      REAL                                                              &
     &  U_LBC_REAL_TEND(LENRIMA(fld_type_u,halo_type_extended,          &
     &                    rima_type_norm),MODEL_LEVELS)                 &
     &, V_LBC_REAL_TEND(LENRIMA(fld_type_v,halo_type_extended,          &
     &                    rima_type_norm),MODEL_LEVELS)                 &
     &, W_LBC_REAL_TEND(LENRIMA(fld_type_p,halo_type_extended,          &
     &                    rima_type_norm),MODEL_LEVELS)                 &
     &, EXNER_LBC_REAL_TEND(LENRIMA(fld_type_p,halo_type_extended,      &
     &                       rima_type_norm),MODEL_LEVELS+1)

! Physics arrays needed by dynamics
      Real                                                              &
     &  rho_km(1-offx:row_length+offx, 1-offy:rows+offy,                &
     &           0:bl_levels-1)                                         &
     &, cH(1-offx:row_length+offx, 1-offy:rows+offy,                    &
     &       model_levels-1)                                            &
     &, wet_to_dry_n (1-offx:row_length+offx, 1-offy:rows+offy,         &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-offx:row_length+offx, 1-offy:rows+offy,       &
     &          model_levels)

! arrays holding information to be passed between physics
! routines.

      Real                                                              &
     &  ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, micro_tends(row_length, rows, bl_levels, 2)
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)

! Radiation fields 1. SW & common with LW.
      Real                                                              &
     &  photosynth_act_rad(row_length, rows)                            &
                                             ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
     &, rad_hr(row_length, rows, bl_levels, 2)                          &
!                                 BL (LW,SW) rad heating rates
     &, surf_radflux(row_length, rows)                                  &
     &, dolr(row_length,rows)                                           &
!       local field "dolr" is distinguished from "dolr_field" (in atm_fields_mod)
                                   ! TOA - surface upward LW
     &, SW_tile(land_field,ntiles)                                      &
                                   ! Surface net SW on land tiles
     &, cos_zenith_angle(row_length, rows)


! mpp-related Arrays
      INTEGER                                                           &
     & g_row_length(0:nproc-1)                                          &
                               ! Table of number of points on a row
     &,g_rows(0:nproc-1)                                                &
                               ! Table number of rows in theta field
     &,g_i_pe(1-halo_i:global_row_length+halo_i) ! processor on my
!               processor-row holding a given value in i direction

! Workspace defined as allocatable arrays, since they each communicate
! fields between near adjacent calls and only need to use memory for
! a subset of the total routine.
      REAL, DIMENSION (:,:,:), ALLOCATABLE::                            &
     &  exner_prime                                                     &
                        ! solution to helmholtz solver
     &, ozone3D                                                         &
                        ! 3d ozone (expanded from zonal) for radiation
     &, dtheta_dr_term                                                  &
     &, depart_lambda, depart_phi, depart_r_theta, depart_r_w

! Allocatable arrays for use in AC_CTL call
      REAL, DIMENSION (:,:,:), ALLOCATABLE::                            &
     &  work_q, work_qcl, work_qcf

      INTEGER                                                           &
     & ErrorStatus      ! Return code : 0 Normal Exit : >0 Error

      CHARACTER*256                                                     &
     & CMessage         ! Error message if return code >0

! Local dynamic arrays:
      REAL, DIMENSION (:), ALLOCATABLE::                                &
!  STASHworki = STASHwork for section i
     & STASHwork1,STASHwork2,STASHwork3,STASHwork4,STASHwork5           &
     &,STASHwork6,STASHwork8,STASHwork9,STASHwork12,STASHwork13         &
     &,STASHwork14,STASHwork17,STASHwork19,STASHwork26,STASHwork30      &
     &,STASHwork10,STASHwork18

! increment diagnostics:
      REAL, DIMENSION (:,:,:), ALLOCATABLE::                            &
     & u_incr_diagnostic, v_incr_diagnostic, T_incr_diagnostic          &
     &,q_incr_diagnostic ,qcl_incr_diagnostic, qcf_incr_diagnostic      &
     &,cf_incr_diagnostic, cfl_incr_diagnostic, cff_incr_diagnostic     &
     &,w_incr_diagnostic

      REAL, DIMENSION (:,:), ALLOCATABLE::                              &
     & w_local_mask
      REAL STASHwork0_dummy(1) ! STASHwork not defined for section 0,
                               !  but required as dummy argument.

! Local arrays for phys1 and phys2 increments for tracers:
      Real :: super_tracer_phys1(1-halo_i:row_length+halo_i,            &
                                 1-halo_j:rows+halo_j,                  &
                                 model_levels, super_array_size)
      Real :: super_tracer_phys2(row_length,                            &
                                 rows,                                  &
                                 model_levels, super_array_size)
      Real :: tracer_phys1(1-halo_i:row_length+halo_i,                  &
                           1-halo_j:rows+halo_j,                        &
                           tr_levels, tr_vars)
      Real :: tracer_phys2(row_length,                                  &
                           rows,                                        &
                           tr_levels, tr_vars)

      Real :: unscaled_dry_rho(1-offx:row_length+offx,                  &
                               1-offy:rows+offy, model_levels)

! Local dynamic arrays for phys1 and phys2 increments for moisture:
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  q_phys1, qcl_phys1, qcf_phys1, q_phys2, qcl_phys2, qcf_phys2    &
     &, qcf2_phys1, qrain_phys1, qgraup_phys1                           &
     &, qcf2_phys2, qrain_phys2, qgraup_phys2                           &
     &, mix_v_phys1, mix_cl_phys1, mix_cf_phys1                         &
     &, mix_v_phys2, mix_cl_phys2, mix_cf_phys2                         &
     &, mix_cf2_phys1, mix_rain_phys1, mix_graup_phys1                  &
     &, mix_cf2_phys2, mix_rain_phys2, mix_graup_phys2                  &
     &, cf_phys1, cfl_phys1, cff_phys1, cf_phys2, cfl_phys2, cff_phys2

! local dynamic arrays for PC2
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  t_inc_pres, q_inc_pres, qcl_inc_pres, qcf_inc_pres, cf_inc_pres &
     &, cfl_inc_pres, cff_inc_pres                                      &
     &, t_dini, q_dini, qcl_dini, qcf_dini, cf_dini                     &
     &, cfl_dini, cff_dini, rhts, qtts, tlts, ptts

!
! Extra variables needed for cycling
! Vars ending in _phys1 are copies holding the value the original
! variable had after exiting phys1.
! Vars ending in _np1 are tn+1 estimates holding the value the original
! variable had at the end of the last cycle (provided that CycleNo>1).
! obtained from the last
! cycle when CycleNo>1.
!
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  R_u_phys1, R_v_phys1, thetastar_phys1, qstar_phys1              &
     &, qclstar_phys1, qcfstar_phys1, qcf2_star_phys1, qrain_star_phys1 &
     &, qgraup_star_phys1, ti_phys1, cca_phys1, area_cld_frac_phys1     &
     &, bulk_cld_frac_phys1, bulk_cld_liq_phys1, bulk_cld_fr_phys1      &
     &, u_np1, v_np1, w_np1, theta_np1, rho_np1, q_np1, qcl_np1         &
     &, qcf_np1, qcf2_np1, qrain_np1, qgraup_np1

      REAL, DIMENSION (:,:), ALLOCATABLE::                              &
     &  z0msea_phys1, zh_phys1, t_land_ctile_phys1, t_sice_ctile_phys1, &
     &  t_surf_phys1, t_sf_tile_phys1, snow_tile_phys1, dolr_phys1,     &
     &  RHO_LBC_REAL_TEND

      INTEGER, DIMENSION (:,:), ALLOCATABLE :: ccb_phys1, cct_phys1

      REAL GS1(land_field)

! local variable
      REAL                                                              &
     &  tot_dry_mass_final                                              &
                            ! mass at end of energy correction period
     &, tot_energy_final                                                &
                            ! energy at end of energy correction period
     &, tot_moist_final                                                 &
                            ! moist at end of energy correction period
     &, energy_corr_now                                                 &
                            ! instanteous energy correction
     &, increment_factor    ! For calculating value of LBC at next TL

! Monsoon variables
      Real                                                              &
     &  lambda_half_width                                               &
     &, phi_half_width                                                  &
     &, p_max                                                           &
     &, p_top                                                           &
     &, p_bottom                                                        &
     &, p_arbitrary                                                     &
     &, lambda_heat_centre                                              &
     &, phi_heat_centre                                                 &
     &, max_heat_per_day                                                &
     &, Mons_newtonian_timescale

!    Local arrays for when using mixing ratios
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  mix_v, mix_cl, mix_cf                                           &
     &, mix_v_star, mix_cl_star, mix_cf_star                            &
     &, mix_v_inter, mix_cl_inter, mix_cf_inter                         &
     &, mix_cf2, mix_rain, mix_graup                                    &
     &, mix_cf2_star, mix_rain_star, mix_graup_star                     &
     &, mix_cf2_inter, mix_rain_inter, mix_graup_inter                  &
     &, q_store, qcl_store, qcf_store                                   &
     &, qcf2_store, qrain_store, qgraup_store
!
! Additional variables needed for cycling when mixing ratios are used
!
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  mix_v_star_phys1, mix_cl_star_phys1, mix_cf_star_phys1          &
     &, mix_cf2_star_phys1, mix_rain_star_phys1, mix_graup_star_phys1   &
     &, mix_v_np1, mix_cl_np1, mix_cf_np1                               &
     &, mix_cf2_np1, mix_rain_np1, mix_graup_np1

!  Suarez-Held variables now declared in CRUNTIMC
      Real, ALLOCATABLE :: RHCPT(:,:,:)
      Integer rhc_row_length
      Integer rhc_rows
      INTEGER    :: minutos ! LOCAL Store value of timestep in min.
!
!  Dummy variables for SCM Diagnostics,
!  for passing down to atmos_physics1 and 2 routines
      INTEGER, PARAMETER :: nSCMDpkgs = 12
      LOGICAL L_SCMDiags(nSCMDpkgs)
      LOGICAL                                                           &
     & L_flux_bc                 ! T if prescribed surface fluxes to be used

      REAL                                                              &
     &  flux_e(row_length, rows)                                        &
                                 ! Surface latent heat flux (W/m^2)
     &, flux_h(row_length, rows)                                        &
                                 ! Surface sensible heat flux (W/m^2)
     &, z0m_scm(row_length, rows)                                       &
                                 ! SCM specified z0m (m)
     &, z0h_scm(row_length, rows)! SCM specified z0m (m)
      REAL, DIMENSION(:,:,:), ALLOCATABLE ::                            &
     &  q_inc_subs, th_inc_subs                                         &
                                  ! subsidence increments
     &, q_inc_ls, th_inc_ls                                             &
                                  ! large scale increments
     &, u_inc_dmp, q_inc_dmp, th_inc_dmp                                &   
                                  !Damping incs
     &, v_inc_dmp
! Tolerance for CycleNo >1
      REAL                                                              &
     &  GCR_run_tol_abs                                                 &
     &, GCR_run_tol_res

! Subgrid turbulence scheme variables.
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     & visc_m                                                           &
!            ! diffusion coefficient for momentum from subgrid
!            ! turbulence scheme
     &,visc_h                                                           &
!            ! diffusion coefficient for heat and moisture from
!            ! subgrid turbulence scheme
     &,FM_3D                                                            &
!            ! stability function for momentum transport.
!            ! level 1 value is dummy for use in diagnostics
     &,FH_3D                                                            &
!            ! stability function for heat and moisture.
!            ! level 1 value is dummy for use in diagnostics
     &,visc_BL_m                                                        &
                     ! visc_m only on BL levels
     &,RNEUTML                                                          &
                     ! mixing length scale (m) (lambda)
     &,shear                                                            &
                     ! S from subgrid turbulence scheme
     &,BL_COEF_KM                                                       &
                     ! RHOKM from BL scheme
     &,BL_COEF_KH    ! RHOKH from BL scheme

       Real                                                             &
     & max_diff  ! max diffusion coeff for run
! Function & Subroutine calls:
      External                                                          &
     & Swap_Bounds                                                      &
     &,Polar_filter                                                     &
     &,Atmos_Physics1                                                   &
     &,STASH                                                            &
     &,SL_Thermo                                                        &
     &,SL_Full_wind                                                     &
     &,Atmos_Physics2                                                   &
     &,dry_static_adj                                                   &
     &,Diff_Divdamp_Ctl                                                 &
     &,pg_update                                                        &
     &,Polar_filter_incs                                                &
     &,Polar_vector_wind_n                                              &
     &,PE_Helmholtz_eul                                                 &
     &,Flux_Rho                                                         &
     &,Polar_Reset_Mean                                                 &
     &,Calc_P_star                                                      &
     &,Calc_Exner_at_theta                                              &
     &,Calc_P_from_Exner                                                &
     &,BOUNDVAL                                                         &
     &,AC_CTL                                                           &
     &,St_diag1                                                         &
     &,St_diag2                                                         &
     &,St_diag3                                                         &
     &,Ereport                                                          &
     &,bottom_w_Calc                                                    &
     &,TFilt_cntl                                                       &
     &,LS_ACF_Brooks                                                    &
     &,set_arcl_dimensions                                              &
     &,set_arcl_clim                                                    &
     &,exppxi               ! Function to extract ppxref info


      Integer                                                           &
     &     exppxi               ! Function to extract ppxref info

      LOGICAL                                                           &
     &  first_atmstep_call                                              &
     &, L_update_lbcs                                                   &
     &, L_apply_lbcs                                                    &
     &, GCR_zero_guess_it

      Integer :: itemp

! time-stepping weights. Values may be different
! at different cycles (UMUI controlled).
      REAL                                                              &
     & alpha1, alpha2, alpha3, alpha4

      ! Local parameters for mixing ratio physics
      Logical, Parameter ::                                             &
                                ! Mixing ratios for atmos_physics_1 and 
                                ! _2 are defined through namelist
     & l_mr_acctl    = .false.                                          &
                                ! Use mixing ratios for ac_ctl
     &,l_mr_qtbalcld = .false.                                          &
                                ! Use mixing ratios for qt_bal_cld
     &,l_mr_tfiltctl = .false.                                          &
                                ! Use mixing ratios for tfilt_ctl
     &,l_mr_pc2      = .false.  ! Use mixing ratios for PC2 routines

      SAVE first_atmstep_call

      DATA first_atmstep_call /.TRUE./

      Logical :: L_physics_store
!
! Oxidant mass-mixing ratios and concentrations, for use in sulphur
! cycle.
      REAL, DIMENSION(:,:,:), ALLOCATABLE ::                            &
     &  O3_MMR                                                          &
     &, H2O2_MMR                                                        &
     &, OH_conc                                                         &
     &, HO2_conc 

! Array to hold atmospheric loss for tracers - kdcorbin, 05/10
   REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: lossrate

   !end step of experiment
   integer :: endstep
!
!- End of header

! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

!! FLUME 
      tr_size = (tr_levels * theta_off_size) ! size of a single tracer

!glr cable 2 lines ----------------------------------------------------
!      L_CABLE = .false. !.true. ! to call CABLE paralell to MOSES 
!      l_cable = .true.  ! to call CABLE instead of MOSES

!      print *,'atm_step before allocation',istep_cur,land_points,sm_levels,NTILES
!      print *,'atm_step: t_i, mype, landpoints',istep_cur, mype,land_points

!     Temporary CABLE initialisation - MRD
      if ( first_atmstep_call ) then
         allocate(SNOW_COND(land_points,NTILES,3))
         allocate(STHU_TILE(land_points,NTILES,sm_levels))
!!$         allocate(STHF_TILE(land_points,NTILES,sm_levels))
         allocate(T_SURF_TILE(land_points,NTILES))
!pfu MOD
         allocate(L_TILE_PTS(land_points,NTILES))
      ALLOCATE(FTL_TILE_CAB(LAND_POINTS,NTILES)   )
      ALLOCATE(FTL_CAB(LAND_POINTS)               )
      ALLOCATE(LE_TILE_CAB(LAND_POINTS,NTILES)    )
      ALLOCATE(LE_CAB(LAND_POINTS)                )
      ALLOCATE(TSTAR_TILE_CAB(LAND_POINTS,NTILES) )
      ALLOCATE(TSTAR_CAB(LAND_POINTS)             )
      ALLOCATE(SMCL_CAB(LAND_POINTS,SM_LEVELS)    )
      ALLOCATE(TSOIL_CAB(LAND_POINTS,SM_LEVELS)   )
      ALLOCATE(USTAR_CAB(LAND_POINTS)                         )
      ALLOCATE(SURF_HTF_CAB(LAND_POINTS)          )
      ALLOCATE(RADNET_TILE(LAND_POINTS,NTILES)    )
      ALLOCATE(CH_CAB(land_points)    )
      ALLOCATE(CD_CAB(land_points)     )
      ALLOCATE(U_S_CAB(land_points)     )
      ALLOCATE(SW_DOWN(row_length,rows)           )
      ALLOCATE(RTSOIL_TILE(land_points,NTILES)    )
      ALLOCATE(GFLUX_TILE(land_points,NTILES)     )
      ALLOCATE(SGFLUX_TILE(land_points,NTILES)    )
      ALLOCATE(HCONS(land_points)                 )
      ALLOCATE(alb_tile(land_points,ntiles,4)     )
      ALLOCATE(surf_down_sw(row_length,rows,4)    )
      ALLOCATE(TOT_ALB(land_points,ntiles)        )
      ALLOCATE(albsoil(land_points)                )
      ALLOCATE(CD(row_length,rows)                )
      ALLOCATE(CHX(row_length,rows)               )
      ALLOCATE(LAND_ALBEDO_CABLE(row_length,rows,4))
      ALLOCATE(TRANSP_TILE(LAND_POINTS,NTILES))
      ALLOCATE(WB_LAKE(LAND_POINTS,NTILES))
      ALLOCATE(TOT_WBLAKE(LAND_POINTS))
      ALLOCATE(TOT_SUBRUN(LAND_POINTS))
      ALLOCATE(wblake_ratio)
         SNOW_COND = -huge(1.)
         STHU_TILE = -huge(1.)
         T_SURF_TILE = -huge(1.)
         L_TILE_PTS = .false.
         if (l_cable) then
          wblake_ratio = 0.0
          ! need init ?
          !WB_LAKE     = 0.0
          !TOT_WBLAKE  = 0.0
          !TOT_SUBRUN  = 0.0
         end if
      end if

       R_w = 0.


! Set non-initialised variables to one

       gi             = 1
       lbc_size       = 1
       ij_v           = 1
       ij_u           = 1
! set Error code to zero
      ErrorStatus = 0

! set the L_Tracer flag
      L_Tracer = ( tr_vars > 0 .OR. L_soot .OR. L_CO2_interactive .OR.  &
     &         L_sulpc_SO2 .OR. L_murk .OR. L_dust .OR. L_biomass .OR.  &
     &         L_USE_CARIOLLE .OR. L_ocff .OR. tr_ukca > 0 )
! timestep information
      timestep = SECS_PER_STEPim(atmos_im) ! timestep in seconds
      radiation_timestep = timestep * A_SW_RADSTEP
      radiation_tstep_diag = timestep * A_SW_RADSTEP_DIAG
      radiation_tstep_prog = timestep * A_SW_RADSTEP_PROG
      timestep_number = STEPim(atmos_im)   ! no. of steps since basis
      l_rad_step = L_SW_RADIATE .OR. L_LW_RADIATE ! Activate radiation
! Activate radiation
      l_rad_step_diag = L_SW_RADIATE_DIAG .OR. L_LW_RADIATE_DIAG
      l_rad_step_prog = L_SW_RADIATE_PROG .OR. L_LW_RADIATE_PROG

      pos_timestep =  timestep
      neg_timestep = -timestep

      IF (L_Backwards) timestep = neg_timestep

      IF (PrintStatus  >=  PrStatus_Normal) THEN

        IF (timestep_number  ==  1 .AND. L_Backwards) THEN
          WRITE(6,*) ''
          WRITE(6,*) '          *************************'
          WRITE(6,*) 'Atm_Step: * INTEGRATING BACKWARDS *'
          WRITE(6,*) '          *************************'
          WRITE(6,*) ''
        ENDIF

        IF (timestep_number  ==  1 .AND. L_RHCPT) THEN
          WRITE(6,*) ''
          WRITE(6,*) 'Atm_Step: Running with diagnostic RHcrit'         &
     &                      //' option active'
          WRITE(6,*) ''
        ENDIF

        IF( L_print_pe .or. mype ==0 ) then
          WRITE(6,*) 'Atm_Step: Timestep ', timestep_number
        ENDIF ! L_print_pe .or. mype ==0

      ENDIF     ! PrintStatus

      IF ((model_domain == mt_global) .and. (L_trivial_trigs)) THEN
        ErrorStatus=123
! DEPENDS ON: ereport
        Call Ereport("ATM_STEP", ErrorStatus,                           &
     &        "Unable to run global model with L_trivial_trigs.")
      ENDIF

      ! Lest 16/12/13 - moved up to allocate loop
      !if (l_cable) then
      ! if (first_atmstep_call) wblake_ratio = 0.0
      !end if

! Set radiation switches in rad_switches_mod
      lrad_ctile_fix    = l_ctile_fix
      lrad_cldtop_t_fix = l_cldtop_t_fix
      lrad_quad_src_fix = l_quad_src_fix
      lrad_ccrad        = l_ccrad
      lrad_3d_ccw       = l_3d_ccw
      lrad_ovrlap       = l_ovrlap
      lrad_ccw_scav     = l_ccw_scav
      lrad_emis_land_gen= l_emis_land_gen
      rad_emis_land_gen = emis_land_gen
! sza
! set cloud inhomegenous option
      if ( .not. L_INHOM_CLOUD ) then
         LRAD_TRIPLECLOUDS = .TRUE.
      endif
! sza end

! Set convection switches in cv_cntl_mod
      lcv_phase_lim     = l_phase_lim
      lcv_3d_cca        = l_3d_cca
      lcv_3d_ccw        = l_3d_ccw
      lcv_ccrad         = l_ccrad


      lcv_pc2_diag_sh   = l_pc2_diag_sh
! set up logical switch for diagnostic printing of l2norms
      L_print_L2norms = .false.
      L_print_L2helm = .false.
      itemp = timestep_number - first_norm_print
      If( print_step > 0 .and. itemp >=0) then
        If( mod( itemp , print_step ) == 0 ) then
          If ( L_diag_L2norms ) then
            if( itemp / print_step  > 11 ) then
         write(6,*)'l2norms printing too often, limited to 12 occasions'
            else
            L_print_L2norms = .true.
            endif  ! itemp / print_step > 11
          endif ! L_diag_L2norms
          If(  L_diag_L2helm ) then
            if( itemp / print_step  > 11 ) then
         write(6,*)'l2norms printing too often, limited to 12 occasions'
            else
            L_print_L2helm = .true.
            endif  ! itemp / print_step > 11
          endif !  L_diag_L2helm
        endif !  mod( itemp , print_step ) == 0)
      endIf ! print_step > 0

! grid information
      If (L_initialise_data .AND.                                       &
     &    first_constant_r_rho_level_new /= -1) Then
        first_constant_r_rho_level = first_constant_r_rho_level_new
      else   !   /=  L_initialise_data
        first_constant_r_rho_level = a_inthd(ih_1_c_rho_level)
      endif  !  L_initialise_data

! Set mpp arrays from parvars.h information

      DO i= 0,nproc-1
        g_row_length(i) = g_lasize(1,fld_type_p,halo_type_no_halo,i)
        g_rows      (i) = g_lasize(2,fld_type_p,halo_type_no_halo,i)
      ENDDO ! i processors

      DO i= 1-halo_i, 0
        g_i_pe(i) = 0
      ENDDO ! i

      DO i= 1,global_row_length
        g_i_pe(i) = g_pe_index_EW(i)
      ENDDO ! i row points

      DO i= global_row_length+1, global_row_length+halo_i
        g_i_pe(i) = nproc_x-1
      ENDDO ! i

! End of Set mpp arrays from parvars.h information

      If (first_constant_r_rho_level  ==  1 ) Then
        first_constant_r_rho_level_m1 = first_constant_r_rho_level
      Else
        first_constant_r_rho_level_m1 = first_constant_r_rho_level - 1
      End If

      lambda_start = datastart(1) - halo_i
      If ( L_regular ) Then
        Do i = 1, row_length
          gi = datastart(1) + i - 1
          lambda_a(i) = (gi - .5) * delta_lambda
        endDo
      else  ! variable resolution
        Do i = 1, row_length
          gi = datastart(1) + i - 1
          lambda_a(i) = glambda_u(gi) - base_lambda
        endDo
      endIf ! L_regular

! set number of levels to check for trajectory inside orography at
! bottom of model, the same parameter is also used inside the
! interpolation routine.

      check_bottom_levels = max(bl_levels,                              &
     &                          interp_vertical_search_tol )

! Define number of work arrays needed in sl_full_wind routine depending
! on model domain choice.
      If (model_domain  ==  mt_global .or.                              &
     &    model_domain  ==  mt_cyclic_lam  .or.                         &
     &    model_domain  ==  mt_bi_cyclic_lam ) Then
        n_Y_arrays    = 1
        n_Yw_arrays   = 1
        n_Yd_arrays   = 1
        n_Ydw_arrays  = 1
      Else If (model_domain  ==  mt_lam) Then
        n_Y_arrays    = 3
        n_Yw_arrays   = 2
        n_Yd_arrays   = 3
        n_Ydw_arrays  = 2
      End If

!  i_start, i_stop, j_start, j_stop define the solution domain
!  j_begin, j_end exclude the poles
      i_start = 1
      i_stop = row_length
      j_start = 1
      j_stop = rows
      j_begin = 1
      j_end = rows
!  mt_bi_cyclic_lam keeps defaults; other domains change as below
      If (model_domain == mt_global) Then
        If (at_extremity(PSouth)) j_begin = 2
        If (at_extremity(PNorth)) j_end = rows - 1
      EndIf  ! model_domain  ==  mt_global
      If (model_domain == mt_lam) Then
        If (at_extremity(PSouth)) then
          j_start = 2
        EndIf ! at_extremity(PSouth)
        If (at_extremity(PNorth)) then
          j_stop = rows - 1
        EndIf ! at_extremity(PNorth)        
        If(at_extremity(PWest)) i_start = 2
        If(at_extremity(PEast)) i_stop = row_length - 1
      End If  ! model_domain  ==  mt_LAM
      If (model_domain == mt_cyclic_lam) Then
        If (at_extremity(PSouth)) then
          j_start = 2
          j_begin = j_start
        EndIf ! at_extremity(PSouth)
        If (at_extremity(PNorth)) then
          j_stop = rows - 1
          j_end = j_stop
        EndIf ! at_extremity(PNorth)        
      End If  ! model_domain ==  mt_cyclic_lam

!  Set variables for surface forcing to false/0.0 for standard UM
      L_flux_bc     = .FALSE.
      flux_e(:,:)  = 0.0
      flux_h(:,:)  = 0.0
! Set variables for roughness length to false/0.0 if not initialising      
      If (.NOT. L_initialise_data) Then
        L_spec_z0     = .FALSE.
        z0m_scm(:,:) = 0.0
        z0h_scm(:,:) = 0.0
      End If
!  Set dummy variables for SCM Diagnostics to false for full UM
      L_SCMDiags(1:nSCMDpkgs)  = .FALSE.
!
! ---------------------------------------------------------------------
! Section 0.1  Initialisation for idealised test problems
!              For standard runs go to end section 0.1
! ---------------------------------------------------------------------

      h_print=0.0

      If (L_subfilter_vert) Then

        ALLOCATE (BL_COEF_KM(1:row_length, 1:rows, bl_levels-1))
        ALLOCATE (BL_COEF_KH(1:row_length, 1:rows, bl_levels-1))

        Do k = 1, BL_LEVELS-1
          Do j = 1, rows
            Do i = 1, row_length
              BL_COEF_KM(i,j,k) = 0.0
              BL_COEF_KH(i,j,k) = 0.0
            End Do
          End Do
        End Do

      Else

        ALLOCATE (BL_COEF_KM(1,1,1))
        ALLOCATE (BL_COEF_KH(1,1,1))

      Endif
      if(L_initialise_data)then

! DEPENDS ON: idl_ni_init
        Call IDL_NI_Init(                                               &
     & R, g, kappa, epsilon, Cp, p_zero, Earth_radius, Pi, two_omega    &
     &,model_domain, row_length, rows, n_rows, model_levels,wet_levels  &
     &,TR_VARS, TR_LEVELS, bl_levels, first_constant_r_rho_level        &
     &,cos_theta_latitude, sec_theta_latitude, f3_at_u, f3_at_v         &
     &,timestep, first_atmstep_call, L_regular                          &
     &,delta_x, delta_y, delta_lambda, delta_phi, base_phi, base_lambda &
     &,A_REALHD(rh_rotlat), A_REALHD(rh_rotlong)                        &
     &, glambda_p(lambda_start), phi_p, glambda_u(lambda_start)         &
     &, phi_v, lambda_p_end, phi_p_end                                  &
     &,r_theta_levels, r_rho_levels, r_at_u, r_at_v, z_orog_print       &
     &,eta_theta_levels, eta_rho_levels                                 &
! Multi-processor
     &,offx,offy,halo_i,halo_j, mype, nproc, at_extremity, datastart    &
     &,gc_all_proc_group, global_row_length, global_rows                &
     &,g_rows, g_row_length, g_datastart, nproc_x                       &
! Primary fields
     &,THETA,RHO, EXNER_THETA_LEVELS                                    &
     &,EXNER_RHO_LEVELS,P                                               &
     &,P_THETA_LEVELS, PSTAR                                            &
     &,Q, QCL, QCF,QCF2,QRAIN                                           &
     &,QGRAUP,CF_BULK, CF_LIQUID                                        &
     &,CF_FROZEN,U, V, W                                                &
     &,U_ADV, V_ADV, W_ADV                                              &
! Lateral boundaries
     &,RIMWIDTHA(rima_type_norm), RIMWEIGHTSA                           &
     &,LENRIMA(1,1,rima_type_norm), LBC_SIZEA(1,1,1,rima_type_norm)     &
     &,LBC_STARTA(1,1,1,rima_type_norm)                                 &
     &,THETA_LBC, THETA_LBC_TEND, EXNER_LBC                             &
     &,EXNER_LBC_TEND, RHO_LBC, RHO_LBC_TEND                            &
     &,Q_LBC, Q_LBC_TEND, QCL_LBC, QCL_LBC_TEND                         &
     &,QCF_LBC, QCF_LBC_TEND, QCF2_LBC                                  &
     &,QCF2_LBC_TEND, QRAIN_LBC, QRAIN_LBC_TEND                         &
     &,QGRAUP_LBC, QGRAUP_LBC_TEND, CF_BULK_LBC                         &
     &,CF_BULK_LBC_TEND, CF_LIQUID_LBC                                  &
     &,CF_LIQUID_LBC_TEND, CF_FROZEN_LBC                                &
     &,CF_FROZEN_LBC_TEND, U_LBC, U_LBC_TEND                            &
     &,V_LBC, V_LBC_TEND, W_LBC, W_LBC_TEND                             &
     &,U_ADV_LBC, U_ADV_LBC_TEND, V_ADV_LBC                             &
     &,V_ADV_LBC_TEND, W_ADV_LBC, W_ADV_LBC_TEND                        &
! Grid info for idealised
     &,A_REALHD(rh_z_top_theta), height_domain, big_layers              &
     &,transit_layers, mod_layers, surface_type, p_surface              &
! Profile settings
     &,tprofile_number, qprofile_number, uvprofile_number, Brunt_Vaisala&
     &,theta_surface, dtheta_dz1, height_dz1, u_in, v_in, height_u_in   &
     &,ujet_lat, ujet_width, u_ramp_start, u_ramp_end, f_plane, r_plane &
     &,q1, num_profile_data, zprofile_data, tprofile_data, qprofile_data&
     &,                  num_uvprofile_data, z_uvprofile_data           &
     &,                  uprofile_data, vprofile_data                   &
     &,max_model_levels, max_num_profile_data, max_num_force_times      &
     &,tforce_option, qforce_option, uvforce_option, num_tforce_levels  &
     &,num_tforce_times, num_qforce_levels, num_qforce_times            &
     &,num_uvforce_levels, num_uvforce_times, z_tforce_data, tforce_data&
     &,z_qforce_data, qforce_data, z_uvforce_data, uforce_data          &
     &,vforce_data, tforce_data_modlev, qforce_data_modlev              &
     &,uforce_data_modlev, vforce_data_modlev                           &
! Dynamical core settings
     &,SuHe_pole_equ_deltaT, SuHe_static_stab                           &
     &,base_frictional_timescale, frictional_timescale                  &
     &,SuHe_sigma_cutoff, SuHe_level_weight, L_SH_Williamson            &
!  Horizontal function parameters
     &,                  t_horizfn_number, uv_horizfn_number            &
     &, t_horizfn_data, L_perturb_t, perturb_magnitude_t                &
     &, L_perturb_q, perturb_magnitude_q, L_perturb_correlate_tq        &
     &, L_perturb_correlate_vert, L_perturb_correlate_time              &
     &, perturb_type, perturb_height                                    &
!  Profiles for fixed lbcs and sponge zones
     &,                  u_ref, v_ref, theta_ref, exner_ref, rho_ref    &
     &,                  q_ref                                          &
     &,L_fix_orog_hgt_lbc, orog_hgt_lbc                                 &
     &,zprofile_orog, idl_interp_option, hf                             &
!  Options
     &,L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2                     &
     &,L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc     &
     &,L_constant_dz, L_rotating, L_fixed_lbcs, L_polar_wind_zero       &
     &,L_wind_balance, L_rotate_winds, L_pressure_balance, L_physics    &
     &,L_dry, L_sponge                                                  &
     &,                  L_trivial_trigs, L_perturb, L_code_test        &
     &,L_cyclone, L_baroclinic                                          &
     &,h_print, timestep_number, h_o_actual, grow_steps, h_o_per_step   &
     &,h_o, grid_number, grid_flat, first_theta_height                  &
     &,thin_theta_height, big_factor, mag, lambda_fraction, phi_fraction&
     &,half_width_x, half_width_y, plat_size_x, plat_size_y, Witch_power&
     &,idl_max_num_bubbles, idl_bubble_option, idl_bubble_max           &
     &,idl_bubble_height, idl_bubble_xoffset, idl_bubble_yoffset        &
     &,idl_bubble_width,  idl_bubble_depth, L_idl_bubble_saturate       &
     &,                  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts       &
     &,                  nproc_y, gc_proc_row_group,n_cca_lev           &
     &,IdlSurfFluxSeaOption,IdlSurfFluxSeaParams,L_flux_bc,flux_h,flux_e&
     &,L_spec_z0, z0m_scm, z0h_scm, roughlen_z0m, roughlen_z0h          &
     &,i_hour, i_minute, i_second                                       &
     &,                  problem_number, rad_hr                         &
     &,OROGRAPHY, TSTAR_TILE, ntiles, land_field, land_index            &
     &,CUMULUS, NBDSC, NTDSC, CCA, CCB                                  &
     &,CCT, CCLWP, TSTAR, LAND, SW_INCS                                 &
     &,LW_INCS, T1_SD, Q1_SD, ZH, CF_AREA                               &
     &,TI, Z0, NTML, U_SEA, V_SEA, U_0_P                                &
     &,                  V_0_P )

        elseif(Problem_number  /=   standard                            &
     &                      .and. timestep_number  ==  1 ) Then
          if (.not. L_Physics .and. mype  ==  0)then
           print*,'Data from dump being used without initialising. '
           print*,'If source of dump is a full model run with orography'
           print*,'then there is no guarantee that run will work'
           print*,'since physics is OFF'
          endif !.not. L_Physics.and. mype  ==  0

        end if   ! L_initialise_data
! ---------------------------------------------------------------------
! End Section 0.1  End of Initialisation for idealised test problems
! ---------------------------------------------------------------------
! DEPENDS ON: timer
      If (Ltimer) Call timer('AS Swap_Bounds',5)

! Update domain halos for time-dependent fields
! DEPENDS ON: set_halos
      call set_halos(                                                   &
     &               U, V, W,                                           &
     &               U_ADV, V_ADV, W_ADV,                               &
     &               THETA, Q, QCL, QCF,                                &
     &               QCF2, QRAIN, QGRAUP,                               &
     &               CF_BULK, CF_LIQUID,                                &
     &               CF_FROZEN,                                         &
     &               RHO, P, P_THETA_LEVELS,                            &
     &               EXNER_RHO_LEVELS,                                  &
     &               EXNER_THETA_LEVELS,                                &
     &               MURK,                                              &
     &               DUST_DIV1,DUST_DIV2,                               &
     &               DUST_DIV3,DUST_DIV4,                               &
     &               DUST_DIV5,DUST_DIV6,                               &
     &               SO2, SO4_AITKEN,                                   &
     &               SO4_ACCU, SO4_DISS, DMS,                           &
     &               NH3, SOOT_NEW, SOOT_AGD,                           &
     &               SOOT_CLD, BMASS_NEW,                               &
     &               BMASS_AGD, BMASS_CLD,                              &
     &               OCFF_NEW, OCFF_AGD, OCFF_CLD,                      &
     &               CO2, TRACER, tracer_ukca,                          &
     &               row_length, rows, n_rows, model_levels, wet_levels,&
     &               offx, offy, halo_i, halo_j, tr_levels, tr_vars,    &
     &               tr_ukca,                                           &
     &               L_MURK,                                            &
     &               L_DUST,                                            &
     &               L_SULPC_SO2, L_SULPC_NH3, L_SULPC_DMS,             &
     &               l_soot, l_biomass, l_ocff, l_co2_interactive,      &
     &               L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,             &
     &               OZONE_TRACER,                                      &
     &               L_USE_CARIOLLE                                     &
     &               )


! DEPENDS ON: timer
      If (Ltimer) Call timer('AS Swap_Bounds',6)

! ---------------------------------------------------------------------
!   Section 0.3  Update lbcs for LAMs
! ---------------------------------------------------------------------
      If ((ErrorStatus == 0) .and. (model_domain == mt_lam)) Then
      
       If ( L_lbc_new ) THEN
         L_update_lbcs = .false.
       Else If ( first_atmstep_call ) THEN
         L_update_lbcs = .false.
       Else If ( RIM_STEPSA == 0 ) THEN
         L_update_lbcs = .false.
       Else If ( L_Fixed_lbcs ) THEN
         L_update_lbcs = .false.
       Else !  Old lbcs and NOT first_atmstep_call 
         L_update_lbcs = .true.
       end if ! L_lbc_new

       If ( L_update_lbcs ) Then

        If (MOD (bndary_offsetim(atmos_im) + stepim(atmos_im),          &
     &           RIM_STEPSA) /= 1 ) Then

          If (RIM_STEPSA  /=  1) Then

! DEPENDS ON: BOUNDVAL
              Call BOUNDVAL(LENRIMA(1,1,rima_type_norm),                &
     &                      L_mcr_qcf2_lbc, L_mcr_qrain_lbc,            &
     &                      L_mcr_qgraup_lbc, L_pc2_lbc,                &
     &                      L_murk_lbc, L_int_uvw_lbc,                  &
     &                      U_LBC, U_LBC_TEND,                          &
     &                      V_LBC, V_LBC_TEND,                          &
     &                      W_LBC, W_LBC_TEND,                          &
     &                      RHO_LBC, RHO_LBC_TEND,                      &
     &                      THETA_LBC, THETA_LBC_TEND,                  &
     &                      Q_LBC, Q_LBC_TEND,                          &
     &                      QCL_LBC, QCL_LBC_TEND,                      &
     &                      QCF_LBC, QCF_LBC_TEND,                      &
     &                      QCF2_LBC, QCF2_LBC_TEND,                    &
     &                      QRAIN_LBC, QRAIN_LBC_TEND,                  &
     &                      QGRAUP_LBC, QGRAUP_LBC_TEND,                &
     &                      CF_BULK_LBC, CF_BULK_LBC_TEND,              &
     &                      CF_LIQUID_LBC, CF_LIQUID_LBC_TEND,          &
     &                      CF_FROZEN_LBC, CF_FROZEN_LBC_TEND,          &
     &                      EXNER_LBC, EXNER_LBC_TEND,                  &
     &                      U_ADV_LBC, U_ADV_LBC_TEND,                  &
     &                      V_ADV_LBC, V_ADV_LBC_TEND,                  &
     &                      W_ADV_LBC, W_ADV_LBC_TEND,                  &
     &                      MURK_LBC, MURK_LBC_TEND,                    &
     &                      TRACER_LBC, TRACER_LBC_TEND,                &
     &                      0, 1, ErrorStatus, CMESSAGE)

           End If       ! RIM_STEPSA  /=  1
         End If ! MOD(bndary_offsetim+stepim,RIM_STEPSA) /= 1 )

        End If !  L_update_lbcs 


        !--------------------------------------------------------------
        ! Idealised UM LBC forcing
        !  If active, update lateral boundary arrays to contain
        !  idealised namelist profile data interpolated in time.
        !--------------------------------------------------------------
        If (L_initialise_data .and. L_force_lbc) Then

! DEPENDS ON: idl_force_lbc
          Call IDL_Force_LBC (                                          &
     &             R, g, Cp, kappa, epsilon, p_zero                     &
     &,            row_length, rows, offx, offy                         &
     &,            halo_i, halo_j, Earth_radius                         &
     &,            LENRIMA(1,1,rima_type_norm)                          &
     &,            timestep, timestep_number                            &
     &,            model_levels, wet_levels                             &
     &,            max_model_levels, max_num_force_times                &
     &,            U_LBC, V_LBC                                         &
     &,            THETA_LBC,Q_LBC                                      &
     &,            U_ADV_LBC,V_ADV_LBC                                  &
     &,            EXNER_LBC                                            &
     &,            r_theta_levels, r_rho_levels                         &
     &,            eta_theta_levels, eta_rho_levels                     &
     &,            height_domain, theta_surface                         &
     &,            pforce_option                                        &
     &,            tforce_option, qforce_option, uvforce_option         &
     &,            num_pforce_times                                     &
     &,            num_tforce_times, num_qforce_times                   &
     &,            num_uvforce_times                                    &
     &,            pforce_time_interval                                 &
     &,            tforce_time_interval, qforce_time_interval           &
     &,            uvforce_time_interval                                &
     &,            p_surface_data                                       &
     &,            tforce_data_modlev, qforce_data_modlev               &
     &,            uforce_data_modlev, vforce_data_modlev               &
     &,            newtonian_timescale )

        End If ! on (L_initialise_data .and. L_force_lbc)

       If ( first_atmstep_call ) THEN
         L_apply_lbcs = .true.
       Else If ( L_lbc_new ) THEN
         L_apply_lbcs = .false.
       Else !  Old lbcs  
         L_apply_lbcs = .true.
       end if ! L_lbc_new

       If ( L_apply_lbcs ) THEN

        !--------------------------------------------------------------
        !           Update primary fields with LAM LBC data
        !--------------------------------------------------------------

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UPDATE_LAM_LBCS',3)

! DEPENDS ON: UPDATE_LAM_LBCS
        CALL UPDATE_LAM_LBCS(                                           &
     &    r_rho_levels, r_theta_levels,                                 &
     &    ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,               &
     &    OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY,                         &
     &    L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,                 &
     &    L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc, &
     &    L_murk, L_murk_lbc,                                           &
     &    L_LBC_balance, L_int_uvw_lbc,                                 &
     &    RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                        &
     &    LENRIMA(1,1,rima_type_norm),                                  &
     &    LBC_SIZEA(1,1,1,rima_type_norm),                              &
     &    LBC_STARTA(1,1,1,rima_type_norm),                             &
     &    THETA_LBC,Q_LBC,QCL_LBC,                                      &
     &    QCF_LBC,QCF2_LBC,QRAIN_LBC,                                   &
     &    QGRAUP_LBC, CF_BULK_LBC,CF_LIQUID_LBC,                        &
     &    CF_FROZEN_LBC, RHO_LBC,EXNER_LBC,                             &
     &    U_LBC,V_LBC,W_LBC,                                            &
     &    U_ADV_LBC,V_ADV_LBC,W_ADV_LBC,                                &
     &    MURK_LBC,                                                     &
     &    THETA,Q,QCL,QCF,                                              &
     &    QCF2,QRAIN,QGRAUP,                                            &
     &    CF_BULK,CF_LIQUID,CF_FROZEN,                                  &
     &    RHO,EXNER_RHO_LEVELS,                                         &
     &    U,V,W,                                                        &
     &    U_ADV,V_ADV,W_ADV,                                            &
     &    MURK,                                                         &
     &    DELTA_PHI, DELTA_LAMBDA,                                      &
     &    BASE_PHI, BASE_LAMBDA,                                        &
     &    TWO_OMEGA, DATASTART,                                         &
     &    lat_rot_NP,                                                   &
     &    GLOBAL_ROW_LENGTH, GLOBAL_ROWS                                &
     &     )

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UPDATE_LAM_LBCS',4)
!
! Must now re-calculate the pressure-based variables, namely pressure
! on both rho and theta levels, exner on theta levels and pstar so that
! they are all consistent with the new LBC-updated values of exner on
! rho levels.
!
! DEPENDS ON: consistent_pressure
        call     Consistent_Pressure (                                  &
     &           exner_rho_levels,                                      &
     &           offx,offy,halo_i,halo_J,                               &
     &           row_length,rows,model_levels,                          &
     &           kappa, g, r_theta_levels, r_rho_levels, rho, p_zero,   &
     &           p, pstar, p_theta_levels,exner_theta_levels)
! ----------------------------------------------------------------------
! Now check that cloud is consistent with moisture fields
! Only really needed in lateral boundary zone but done everywhere
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
!  Check liquid/ice cloud fraction is zero if qcl/qcf = 0
! ----------------------------------------------------------------------
!

! DEPENDS ON: cloud_check
        Call Cloud_check(                                               &
     &                 j_start, j_stop, i_start, i_stop                 &
     &,                rows, row_length, wet_levels                     &
     &,                halo_i, halo_j                                   &
     &,                QCL, QCF                                         &
     &,                CF_LIQUID, CF_FROZEN                             &
     &,                CF_AREA, CF_BULK)

        End If ! L_apply_lbcs

        if ( L_int_uvw_lbc ) then

! Obtain tendencies in the boundary zone for u, v, w and use to make.
! lbcs for u_adv, v_adv, w_adv,
! This means that  u_adv, v_adv and w_adv
! can be removed from the lateral boundary files.
! Similar code is needed later for tendencies for the solver

          L_do_halos=.TRUE.
          L_do_boundaries=.TRUE.

          If (RIM_STEPSA  ==  0) Then
            increment_factor=0.0
          Else
            increment_factor=1.0/                                       &
     &      (RIM_STEPSA-MOD(Timestep_Number-1,RIM_STEPSA))
          End If

          lbc_size=LENRIMA(fld_type_p,halo_type_extended,               &
     &             rima_type_norm)

          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              W_LBC_REAL_TEND(i,k) = W_LBC(i,k) +                       &
     &                                   0.5 * increment_factor *       &
     &                           ( W_LBC_TEND(i,k) - W_LBC(i,k) )
            END DO
          END DO

! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH, ROWS, HALO_I, HALO_J,                             &
     &    MODEL_LEVELS, fld_type_p, W_ADV,                              &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    W_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm), RIMWIDTHA(rima_type_norm),         &
     &    RIMWEIGHTSA, AT_EXTREMITY,                                    &
     &    L_do_boundaries, L_do_halos)

          lbc_size=LENRIMA(fld_type_u,halo_type_extended,               &
     &             rima_type_norm)
 
          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              U_LBC_REAL_TEND(i,k) = U_LBC(i,k) +                       &
     &                                      0.5 * increment_factor *    &
     &                          ( U_LBC_TEND(i,k) - U_LBC(i,k) )
            END DO
          END DO

! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH, ROWS, HALO_I, HALO_J,                             &
     &    MODEL_LEVELS, fld_type_u, U_ADV,                              &
     &    LENRIMA(fld_type_u,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_u,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_u,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    U_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm), RIMWIDTHA(rima_type_norm),         &
     &    RIMWEIGHTSA, AT_EXTREMITY,                                    &
     &    L_do_boundaries, L_do_halos)
     
          lbc_size=LENRIMA(fld_type_v,halo_type_extended,               &
     &             rima_type_norm)
 
          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              V_LBC_REAL_TEND(i,k) = V_LBC(i,k) +                       &
     &                                    0.5 * increment_factor *      &
     &                            ( V_LBC_TEND(i,k) - V_LBC(i,k) )
            END DO
          END DO

! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH, N_ROWS, HALO_I, HALO_J,                           &
     &    MODEL_LEVELS, fld_type_v, V_ADV,                              &
     &    LENRIMA(fld_type_v,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_v,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_v,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    V_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm), RIMWIDTHA(rima_type_norm),         &
     &    RIMWEIGHTSA, AT_EXTREMITY,                                    &
     &    L_do_boundaries, L_do_halos)

        endif !  L_int_uvw_lbc
          

      ENDIF     !   model_domain  ==  mt_lam
!QAN - No code for mt_cyclic_lam yet
       
! ----------------------------------------------------------------------
! Section 0.1  Filter winds and theta near poles if active
!              Do horizontal diffusion as a filter if active
! ----------------------------------------------------------------------
!    diagnostic printing of l2norms
      If( L_print_L2norms ) then
        If( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms start of Timestep ***'
        End If ! L_print_pe .or. mype ==0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .true., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .true., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      End If !  L_print_L2norms

! Call timer for diffusion code
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('NI_filter_Ctl',3)

      If ( L_filter ) then

! section 13:
        IF( SF(0,13) ) THEN    ! Diagnostics required for this section

! Allocate diagnostic space for STASH
          ALLOCATE (STASHwork13(STASH_maxlen(13,A_im)))
        
        ENDIF

! DEPENDS ON: ni_filter_ctl
       Call NI_filter_Ctl(                                              &
     &                      THETA,                                      &
     &                      U, V, W,                                    &
     &                      EXNER_RHO_LEVELS, RHO,                      &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      r_theta_levels, r_rho_levels,               &
     &                      r_at_u, r_at_v, delta_lambda, delta_phi,    &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      cos_v_latitude, sec_v_latitude,             &
     &                      polar_filter_north_lat_limit,               &
     &                      polar_filter_south_lat_limit,               &
     &                      polar_filter_coefficient,                   &
     &                      polar_filter_n_sweeps,                      &
     &                      polar_filter_step_per_sweep,                &
     &                      polar_filter_lat_limit,                     &
     &                      max_121_rows, u_sweeps, v_sweeps,           &
     &                      global_u_filter, global_v_filter,           &
     &                      u_begin, u_end, v_begin, v_end,             &
     &                      diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
     &                      diff_coeff_thermo, diff_coeff_wind,         &
     &                      diff_order_thermo, diff_order_wind,         &
     &                      first_constant_r_rho_level,                 &
     &                      first_constant_r_rho_level_m1,              &
     &                      top_filt_start, top_filt_end,               &
     &                      up_diff, max_updiff_levels,                 &
     &                      horizontal_level, mype,                     &
     &                      global_row_length, global_rows,             &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      nproc, nproc_x, nproc_y, datastart,         &
     &                      neighbour, at_extremity, model_domain,      &
     &                      gc_proc_row_group, gc_proc_col_group,       &
     &                      L_polar_filter, L_pofil_new,                &
     &                      L_pfcomb, L_pftheta, L_pfuv,                &
     &                      L_pfw, L_pfexner, L_diff_exner,             &
     &                      L_diff_thermo, L_diff_wind, L_diff_w,       &
     &                      L_pofil_hadgem2, Ltimer,exner_theta_levels, &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &                      STASHwork13)

      if(L_pfexner .and. L_pofil_new)then
! DEPENDS ON: consistent_pressure
        call     Consistent_Pressure (                                  &
     &           exner_rho_levels,                                      &
     &           offx,offy,halo_i,halo_J,                               &
     &           row_length,rows,model_levels,                          &
     &           kappa, g, r_theta_levels, r_rho_levels, rho, p_zero,   &
     &           p, pstar, p_theta_levels,exner_theta_levels)
      endif  !  (L_pfexner .and. l_pofil_new)

      End If    !  L_filter

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('NI_filter_Ctl',4)

! ----------------------------------------------------------------------
! Section 0.2  If free-slip then put non-zero w field on bottom boundary
! ----------------------------------------------------------------------
      if ( L_free_slip ) then
! DEPENDS ON: bottom_w_calc
        Call Bottom_w_Calc(                                             &
     &                     r_theta_levels,                              &
     &                     u, v,                                        &
     &                     w, w_adv,                                    &
     &                     sec_theta_latitude, delta_lambda, delta_phi, &
     &                     glambda_p(lambda_start), phi_p,              &
     &                     glambda_u(lambda_start), phi_v,              &
     &                     rows, n_rows, row_length, model_levels,      &
     &                     model_domain, L_regular,                     &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     gc_proc_row_group, at_extremity,             &
     &                     global_row_length)

! DEPENDS ON: swap_bounds
        call Swap_Bounds                                                &
     &                  (W,                                             &
     &                   row_length, rows, model_levels+1,              &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
        call Swap_Bounds                                                &
     &                  (W_ADV,                                         &
     &                   row_length, rows, model_levels+1,              &
     &                   halo_i, halo_j, fld_type_p, .false.)

      endif !( L_free_slip )
! Call to the RANDOM PARAMETERS (STPH_RP) subroutine
       IF (L_Physics .and. L_RP) THEN
        minutos=timestep/60
        IF(mod(i_hour,3) == 0 .and. i_minute == minutos) THEN
         IF (PrintStatus  >=   PrStatus_Normal) THEN
            WRITE(6,*) 'CALLING RANDOM PARAMETERS'
         ENDIF
! DEPENDS ON: stph_rp
         CALL stph_rp(max_model_levels,                                 &
     &            RHCRIT,RHCRIT_max,RHCRIT_min,                         &
     &            CI_max,CI_min,                                        &
     &            GWD_FRC,GWD_FRC_max,GWD_FRC_min,                      &
     &            KAY_GWAVE,KAY_GWAVE_max,KAY_GWAVE_min,                &
     &            par_mezcla_max,par_mezcla_min,                        &
     &            G0_max,G0_min,                                        &
     &            G0_RP,par_mezcla)
        ELSE
         IF (PrintStatus  >=  PrStatus_Normal) THEN
             WRITE(6,*) 'NOT CALLING RANDOM PARAMETERS'
             WRITE(6,*) 'This routine is only called every 3hrs'
         ENDIF
        ENDIF
      ENDIF
! Call to the RANDOM PARAMETERS2 (STPH_RP2) subroutine
       IF (L_Physics .and. L_RP2) THEN
        minutos=timestep/60
        IF(mod(i_hour,3) == 0 .and. i_minute == minutos) THEN
         IF (PrintStatus  >=   PrStatus_Normal) THEN
            WRITE(6,*) 'CALLING RANDOM PARAMETERS2'
         ENDIF
! DEPENDS ON: stph_rp2
          CALL stph_rp2(max_model_levels,                               &
     &            RHCRIT,RHCRIT_max,RHCRIT_min,                         &
     &            CI_max,CI_min,                                        &
     &            GWD_FRC,GWD_FRC_max,GWD_FRC_min,                      &
     &            KAY_GWAVE,KAY_GWAVE_max,KAY_GWAVE_min,                &
     &            par_mezcla_max,par_mezcla_min,                        &
     &            G0_max,G0_min,G0_RP,par_mezcla,M_CI,                  &
     &            M_CI_max,M_CI_min,                                    &
     &            Charnock,Charnock_max,Charnock_min)
        ELSE
         IF (PrintStatus  >=  PrStatus_Normal) THEN
             WRITE(6,*) 'NOT CALLING RANDOM PARAMETERS2'
             WRITE(6,*) 'This routine is only called every 3hrs'
         ENDIF
        ENDIF
      ENDIF
! ----------------------------------------------------------------------
!    diagnostic printing of l2norms
      If( L_print_L2norms ) then
        If( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after Filter_ctl ***'
        End If ! L_print_pe .or. mype ==0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .true., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .true., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      End If !  L_print_L2norms
! ---------------------------------------------------------------

! ---------------------------------------------------------------
!    diagnostic printing at timestep 1
! ---------------------------------------------------------------
       if( L_diag_print .and. timestep_number ==1 ) then
! DEPENDS ON: print_diag
          Call Print_diag(                                              &
     &                  U, V, THETA,                                    &
     &                  RHO, W,                                         &
     &                  Q, QCL, QCF,                                    &
     &                  rows, n_rows, row_length,                       &
     &                  model_levels, wet_levels, model_domain,         &
     &                  global_row_length, global_rows,                 &
     &                  r_theta_levels, r_rho_levels,                   &
     &                  r_at_u, r_at_v , Pi,                            &
     &                  FV_sec_theta_latitude, FV_cos_theta_latitude,   &
     &                  cos_v_latitude,                                 &
     &                  cos_theta_longitude, sin_theta_longitude,       &
     &                  offx, offy, halo_i, halo_j,                     &
     &                  mype, nproc, at_extremity, datastart,           &
     &                  gc_proc_row_group, delta_lambda, delta_phi,     &
     &                  timestep_number, 1, diag_interval,              &
     &                  rpemax, rpemin, ipesum, rpesum, w_print_limit,  &
     &                  L_print_pe, L_print_w,                          &
     &                  L_print_wmax, L_print_max_wind,                 &
     &                  L_print_div, L_print_lapse, L_print_theta1,     &
     &                  L_print_shear, L_diag_wind, L_diag_noise,       &
     &                  max_w_run, max_wind_run, min_theta1_run,        &
     &                  dtheta1_run, max_div_run, min_div_run,          &
     &                  min_lapse_run, max_shear_run, time_max_shear,   &
     &                  time_div_max, time_div_min, time_lapse_min,     &
     &                  time_w_max, time_max_wind, time_theta1_min,     &
     &                  max_KE_run, min_KE_run, max_noise_run,          &
     &                  time_KE_max, time_KE_min, time_noise_max )
! DEPENDS ON: um_fort_flush
         if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
       endif     !  timestep_number ==1
! ---------------------------------------------------------------
! Section 1.0  Call Atmospheric Physics1
! ----------------------------------------------------------------------
      
      !
      ! Biogenic aerosol climatology for the climate and NWP models
      !
      IF (L_USE_BIOGENIC) THEN 
        do k=1, model_levels
          do j=1, rows
            do i=1, row_length
              biogenic(i,j,k) = ARCLBIOG_BG(i,j,k)
            end do
          end do
        end do
      END IF
            
      !
      ! Aerosol climatologies for the Numerical Weather Prediction
      ! model. Model switches and climatologies are gathered into 
      ! bigger arrays.
      !
      ! First, set the internal model switches according to the
      ! value of the CNTLATM switches, and determine how many 
      ! components we need.
      !
!DEPENDS ON: set_arcl_dimensions
      CALL set_arcl_dimensions(                                         &
     &                         L_USE_ARCLBIOM,                          &
     &                         L_USE_ARCLBLCK,                          &
     &                         L_USE_ARCLSSLT,                          &
     &                         L_USE_ARCLSULP,                          &
     &                         L_USE_ARCLDUST,                          &
     &                         L_USE_ARCLOCFF,                          &
     &                         L_USE_ARCLDLTA,                          &
     &                         n_arcl_species,                          &
     &                         n_arcl_compnts,                          &
     &                         L_USE_ARCL                               &
     &                        )
      
      !
      ! If the aerosol climatology for NWP is used, n_arcl_species
      ! is larger than 0. In that case, allocate the array gathering
      ! component mass-mixing ratio and take the values from the
      ! arrays in arg_atm_fields.h.
      !
      if (n_arcl_species > 0) Then
      
        allocate(arcl(row_length, rows, model_levels, n_arcl_compnts))
      
!DEPENDS ON: set_arcl_clim
        call set_arcl_clim(                                             &
                           ! Array dimensions
     &                     row_length, rows, model_levels,              &
     &                     n_arcl_compnts,                              &
                           ! Internal model switches
     &                     L_USE_ARCL,                                  &
                           ! Climatologies from ancillary files
     &                     arclbiom_fr, arclbiom_ag, arclbiom_ic,       &
     &                     arclblck_fr, arclblck_ag,                    &
     &                     arclsslt_fi, arclsslt_jt,                    &
     &                     arclsulp_ac, arclsulp_ak, arclsulp_di,       &
     &                     arcldust_b1, arcldust_b2, arcldust_b3,       &
     &                     arcldust_b4, arcldust_b5, arcldust_b6,       &
     &                     arclocff_fr, arclocff_ag, arclocff_ic,       &
     &                     arcldlta_dl,                                 &
                           ! Internal climatology array
     &                     arcl,                                        &
                           ! Component array indices
     &                     i_arcl_compnts                               &
     &       )
        
      else
         allocate ( arcl(1,1,1,1) )
      end if
            
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS1 Atmos_Phys1',5)

      IF (L_Physics .AND. ErrorStatus == 0) THEN

! Set logicals for tropopause diagnostics

        L_O3_trop_level  = ( sf(280,2) ) .AND.                          &
     &                     ( ( i_ozone_int  ==  IO3_TROP_MAP) .OR.      &
     &                       ( i_ozone_int  ==  IO3_TROP_MAP_MASSCON) )
        L_O3_trop_height = ( sf(281,2) ) .AND.                          &
     &                     ( ( i_ozone_int  ==  IO3_TROP_MAP) .OR.      &
     &                       ( i_ozone_int  ==  IO3_TROP_MAP_MASSCON) )
        L_T_trop_level   = ( sf(282,2) ) .AND.                          &
     &                     ( ( i_ozone_int  ==  IO3_TROP_MAP) .OR.      &
     &                       ( i_ozone_int  ==  IO3_TROP_MAP_MASSCON) )
        L_T_trop_height  = ( sf(283,2) ) .AND.                          &
     &                     ( ( i_ozone_int  ==  IO3_TROP_MAP) .OR.      &
     &                       ( i_ozone_int  ==  IO3_TROP_MAP_MASSCON) )

! Allocate space for tropopause diagnostics
        IF( L_O3_trop_level ) THEN
          ALLOCATE ( O3_trop_level(row_length,rows) )
        ELSE
          ALLOCATE ( O3_trop_level(1,1) )
        ENDIF
        IF( L_O3_trop_height ) THEN
          ALLOCATE ( O3_trop_height(row_length,rows) )
        ELSE
          ALLOCATE ( O3_trop_height(1,1) )
        ENDIF
        IF( L_T_trop_level ) THEN
          ALLOCATE ( T_trop_level(row_length,rows) )
        ELSE
          ALLOCATE ( T_trop_level(1,1) )
        ENDIF
        IF( L_T_trop_height ) THEN
          ALLOCATE ( T_trop_height(row_length,rows) )
        ELSE
          ALLOCATE ( T_trop_height(1,1) )
        ENDIF

! Check whether ozone option is applicable to the specified
! ozone ancillary file. If not output error message.
        IF (PrintStatus  >=   PrStatus_Diag) THEN   
            WRITE(6,*) 'Atm Step: Lexpand_ozone', Lexpand_ozone
        END IF
  
        IF ((Lexpand_ozone).and.(I_ozone_int == 1)) THEN
           ErrorStatus=123
! DEPENDS ON: ereport
           Call Ereport("ATM_STEP", ErrorStatus,                        &
     &        "A 2D ozone ancillary has been specified with a" //       &
     &        "3D ozone option.")
        ENDIF

        IF (.not.Lexpand_ozone) THEN
           IF (I_ozone_int == 2) THEN
              ErrorStatus=123
! DEPENDS ON: ereport
              Call Ereport("ATM_STEP", ErrorStatus,                     &
     &           "A 3D ozone ancillary has been specified with a" //    &
     &           "2D ozone option.")
           ENDIF

           IF (I_ozone_int == 5) THEN
              ErrorStatus=123
! DEPENDS ON: ereport
              Call Ereport("ATM_STEP", ErrorStatus,                     &
     &           "A 3D ozone ancillary has been specified with a" //    &
     &           "2D ozone option.")
           ENDIF
        ENDIF


!     Convert the ozone field supplied from the dump to a full 3-D
!     array. At 5.2 we require the ozone field to be specified on
!     the same number of levels as the model's grid: this will later
!     be generalized to allow for a different number of levels.
        If (Lexpand_ozone) then
          nd_o3=rows*ozone_levels
        Else
          nd_o3=row_length*rows*ozone_levels
        Endif
!  STOCHEM size
        nd_stochem=row_length*rows*model_levels

! Allocate workspace, just required in atmos_physics1
        ALLOCATE (ozone3D (row_length, rows, ozone_levels) )
! DEPENDS ON: o3_to_3d
        Call O3_to_3D(                                                  &
     &                lexpand_ozone, i_ozone_int,                       &
     &                rows, row_length, model_levels, ozone_levels,     &
     &                halo_i, halo_j, offx, offy, at_extremity,         &
     &                a_realhd(rh_z_top_theta),                         &
     &                THETA,                                            &
     &                r_theta_levels, eta_theta_levels,                 &
     &                EXNER_THETA_LEVELS,                               &
     &                RHO,                                              &
     &                r_rho_levels, eta_rho_levels,                     &
     &                EXNER_RHO_LEVELS,                                 &
     &                nd_o3, O3,                                        &
     &                l_use_stochem_o3,                                 &
     &                nd_stochem,O3_STOCH,                              &
     &  min_trop_level, max_trop_level,                                 &
     &  L_O3_trop_level,L_O3_trop_height,                               &
     &  L_T_trop_level,L_T_trop_height,                                 &
     &  O3_trop_level,O3_trop_height,                                   &
     &  T_trop_level,T_trop_height,                                     &
     &  gc_proc_row_group,                                              &
     &  global_row_length,                                              &
     &                ozone3D,                                          &
     &                ErrorStatus, cmessage                             &
     &                )
!
!     If the ozone_tracer is initialised to 0 then it should be reset to something realistic
!     so set it to climatology expanded to 3D
     
       IF (L_USE_CARIOLLE) THEN
          IF (OZONE_TRACER(1,1,1) == 0.0) THEN
             WRITE(6,*)'O3 tracer must not be set to 0 reset to clim'
        
             DO k=1,model_levels
               DO j=1,rows
                 DO i=1,row_length        
                   OZONE_TRACER(i,j,k)=ozone3D(i,j,k)      
                 END DO
               END DO
             END DO
         
! Halos updated
! DEPENDS ON: swap_bounds
             call Swap_Bounds(OZONE_TRACER,                             &
                        row_length, rows, model_levels,                 &
                        offx, offy, fld_type_p, .false.)
!DEPENDS ON: fill_external_halos
             CALL FILL_EXTERNAL_HALOS(OZONE_TRACER,row_length, rows,    &
                          model_levels,offx,offy)
          
          ELSE  
            WRITE(6,*) 'At least first row of Ozone tracer is not zero'
          END IF ! End of check for zero ozone
       END IF
!     After a successful call ozone3D will contain the mixing ratios of
!     ozone on the model's grid. Uses ozone tracer from the last timestep
!     in the radiation scheme. Radiation is done first so all variables used 
!     are from the last timestep.
!    

       IF (L_USE_OZONEINRAD) THEN                           
             DO k = 1, model_levels                           
               DO j = 1, rows                             
                 DO i = 1, row_length                                  
                   temp_ozone(i,j,k)=OZONE_TRACER(i,j,k)

!  Use the 3D ozone in radiation but check for zero or -ve ozone as some 
!  tropospheric analysed ozone may be -ve; do not use -ve values.
                     IF (temp_ozone(i,j,k) > 0.0) THEN   
                        ozone3D(i,j,k) = temp_ozone(i,j,k)
                     ELSE   
                       IF (j > 1 .and. j < rows .and.     &
     &                    i > 1 .and. i < row_length)       THEN
                         IF (OZONE_TRACER(i,j-1,k) > 0.0)   THEN
                           ozone3D(i,j,k) = OZONE_TRACER(i,j-1,k)
                         ELSE
                           IF (OZONE_TRACER(i,j+1,k) > 0.0) THEN
                             ozone3D(i,j,k) = OZONE_TRACER(i,j+1,k)
                           ELSE
                             IF (OZONE_TRACER(i-1,j,k) > 0.0)    THEN
                              ozone3D(i,j,k) = OZONE_TRACER(i-1,j,k)
                             ELSE
                                IF (OZONE_TRACER(i+1,j,k) > 0.0) THEN
                                   ozone3D(i,j,k) = OZONE_TRACER(i+1,j,k)
                                END IF   ! Check on i+1>0
                             END IF    ! Check on i-1>0
                           END IF    ! Check on j+1>0
                         END IF    ! Check on j-1>0
                       END IF    ! Check on j (1:rows) and i (1:row_length)
                     END IF    ! Check on temp_ozone > 0                                          
                 END DO                                    
               END DO                              
             END DO

       END IF                                                    
      
      END IF !  L_Physics
!
! Set LAND_PTS_TRIF and NPFT_TRIF according to TRIFFID on/off
! and set pointers to single or multi pool soil carbon, and
! soil carbon dimensions to be used in subroutines.
      IF (L_TRIFFID) THEN
        LAND_PTS_TRIF = LAND_FIELD
        NPFT_TRIF = NPFT
        CS => SOIL_CARB1
        RSA => RSP_S_ACC1
        DIM_CS1 = 4
        DIM_CS2 = LAND_FIELD
      ELSE
        LAND_PTS_TRIF = 1
        NPFT_TRIF = 1
        CS => SOIL_CARB
        RSA => RSP_S_ACC
        DIM_CS1 = 1
        DIM_CS2 = 1
      ENDIF
!     dim_cs2 needs to be modified for Carbon fluxes
      IF ( l_cable ) DIM_CS2 = LAND_FIELD
     ! Lestevens 23apr13
      IF ( l_cable ) LAND_PTS_TRIF = LAND_FIELD
      IF ( l_cable ) NPFT_TRIF = NPFT
!
!  set up CO2 field to be passed down
!
      IF (L_CO2_INTERACTIVE) THEN
        CO2_DIM_LEN = row_length
        CO2_DIM_ROW = rows
        CO2_DIM_LEV = model_levels
      ELSE
        CO2_DIM_LEN = 1
        CO2_DIM_ROW = 1
        CO2_DIM_LEV = 1
      ENDIF

      if (L_use_seasalt_autoconv .OR. L_use_seasalt_sulpc .OR.          &
     &         L_use_seasalt_indirect .OR. L_use_seasalt_direct) then
        salt_dim1=row_length
        salt_dim2=rows
        salt_dim3=model_levels
      else
        salt_dim1=1
        salt_dim2=1
        salt_dim3=1
      endif

      if (L_use_seasalt_sulpc .OR. L_DMS_em_inter                       &
     &                        .OR. L_DMS_Ointer) then
        aero_dim1=row_length
        aero_dim2=rows
      else
        aero_dim1=1
        aero_dim2=1
      endif

      if (L_use_seasalt_sulpc) then
        aero_dim3=model_levels
      else
        aero_dim3=1
      endif

! Allocate additional microphysics variables to full size
! if in use, otherwise allocate minimum amount of space

      If (L_mcr_qcf2) Then  ! Allocate second cloud ice
        Allocate ( qcf2_star(1-offx:row_length+offx,                    &
     &                       1-offy:rows+offy, wet_levels) )
      Else
        Allocate ( qcf2_star(1,1,1) )
      End If

      If (L_mcr_qrain) Then  ! Allocate rain
        Allocate ( qrain_star(1-offx:row_length+offx,                   &
     &                        1-offy:rows+offy, wet_levels) )
      Else
        Allocate ( qrain_star(1,1,1) )
      End If

      If (L_mcr_qgraup) Then  ! Allocate graupel
        Allocate ( qgraup_star(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
      Else
        Allocate ( qgraup_star(1,1,1) )
      End If

! The _star fields are used to store the increments to theta, q, qcl,
! and qcf
      IF (L_Physics .AND. ErrorStatus == 0) THEN

! Allocate diagnostic space for STASH
      If (L_radiation) then
        ALLOCATE (STASHwork1(STASH_maxlen(1,A_im)))
        ALLOCATE (STASHwork2(STASH_maxlen(2,A_im)))
      endif
      If (L_rain) then
        ALLOCATE (STASHwork4(STASH_maxlen(4,A_im)))
      endif
      IF(L_gwd .or. l_use_ussp)then
        ALLOCATE (STASHwork6(STASH_maxlen(6,A_im)))
      endif
        IF (L_EMCORR) THEN   ! only if energy correction required
          ALLOCATE (STASHwork14(STASH_maxlen(14,A_im)))
        ENDIF   ! l_emcorr

      ! -----------------------------------------------------------
      ! Convert to mixing ratios from specific humidities if needed
      ! -----------------------------------------------------------
      If (l_mr_physics1) Then

        ! Allocate q_store variables

        allocate ( q_store  (1-halo_i:row_length+halo_i,                &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        allocate ( qcl_store(1-halo_i:row_length+halo_i,                &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        allocate ( qcf_store(1-halo_i:row_length+halo_i,                &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        If (L_mcr_qcf2) then
          allocate ( qcf2_store(1-halo_i:row_length+halo_i,             &
     &                       1-halo_j:rows+halo_j, wet_levels) )
        Else
          allocate ( qcf2_store(1,1,1) )
        End if
        If (L_mcr_qrain) then
          allocate ( qrain_store(1-halo_i:row_length+halo_i,            &
     &                        1-halo_j:rows+halo_j, wet_levels) )
        Else
          allocate ( qrain_store(1,1,1) )
        End if
        If (L_mcr_qgraup) then
          allocate ( qgraup_store(1-halo_i:row_length+halo_i,           &
     &                         1-halo_j:rows+halo_j, wet_levels) )
        Else
          allocate ( qgraup_store(1,1,1) )
        End if

        ! Hold a copy of d1(q) etc. in q_store to place back in d1
        ! after the conversions are finished

        Do k = 1, wet_levels
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              q_store(i,j,k)   = q(i,j,k)   ! Vapour
              qcl_store(i,j,k) = qcl(i,j,k) ! Liquid
              qcf_store(i,j,k) = qcf(i,j,k) ! Ice
            End do
          End do
        End do

        If (L_mcr_qcf2) then  ! Ice2
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qcf2_store(i,j,k) = qcf2(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qrain) then  ! Rain
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qrain_store(i,j,k) = qrain(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qgraup) then  ! Graupel
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qgraup_store(i,j,k) = qgraup(i,j,k)
              End do
            End do
          End do
        End if

        ! Convert d1 values to mixing ratios, including the halos,
        ! using the _store values as inputs (d1s are outputs)

! DEPENDS ON: q_to_mix_halo
        call q_to_mix_halo (row_length, rows, wet_levels,               &
     &           halo_i, halo_j,                                        &
     &           q_store, qcl_store, qcf_store,                         &
     &           qcf2_store, qrain_store, qgraup_store,                 &
     &           L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                 &
     &           Q, QCL, QCF,                                           &
     &           QCF2, QRAIN, QGRAUP                                    &
     &          )

        ! d1 values now contain mixing ratios (including halos).
        ! _store values are specific humidities.

      ! -----------------------------------------------------------
      ! End of conversion to mixing ratios for atmos_physics1 call
      ! -----------------------------------------------------------
      End if  ! l_mr_physics1

        If (L_pc2) Then
!         Calculate the intital Relative humidity wrt TL (rhts)
          allocate(rhts(row_length,rows,wet_levels))
! DEPENDS ON: pc2_rhtl
          Call pc2_rhtl(halo_i, halo_j, offx, offy                      &
     &,     wet_levels,row_length,rows                                  &
     &,     THETA, EXNER_THETA_LEVELS                                   &
     &,     Q,QCL,P_THETA_LEVELS,rhts                                   &
     &,     l_mr_physics1)
! Cusack interpolation needs RH on large_levels from start of timestep
          allocate(tlts(row_length,rows,wet_levels))
          allocate(qtts(row_length,rows,wet_levels))
          allocate(ptts(row_length,rows,wet_levels))
	  do k=1,wet_levels
            do j=1,rows
              do i=1,row_length
                 tlts(i,j,k)=theta(i,j,k)*exner_theta_levels(i,j,k)     &
     &              -QCL(i,j,k)*LC/CP
                 ptts(i,j,k)=p_theta_levels(i,j,k)
                 qtts(i,j,k)=Q(i,j,k)+QCL(i,j,k)
              end do
            end do
          end do
        EndIf

! NB: the star variables and R_u and R_v have not been set in the
!     halo region yet.

! DEPENDS ON: timer
        If (Ltimer) Call timer('Atmos_Physics1',3)

           do l=1,land_field
              frac_control(l,1)=FRAC_CON1(l)
              frac_control(l,2)=FRAC_CON2(l)
              frac_control(l,3)=FRAC_CON3(l)
              frac_control(l,4)=FRAC_CON4(l)
              frac_control(l,5)=FRAC_CON5(l)
              frac_control(l,6)=FRAC_CON6(l)
              frac_control(l,7)=FRAC_CON7(l)
              frac_control(l,8)=FRAC_CON8(l)
              frac_control(l,9)=FRAC_CON9(l)
           enddo
          
      istep_cur = istep_cur + 1  ! For CABLE
! NB if you are changing the argument list to atmos_physics1, please
! do an equivalent change in routine scm_main to keep the single column
! model consistent.

   CALL iday_kick(i_day_number, cable_lai)

! DEPENDS ON: atmos_physics1
        Call Atmos_Physics1(                                            &
! Parallel variables
     & halo_i, halo_j, offx, offy, global_row_length, global_rows       &
     &,gc_proc_row_group,gc_proc_col_group, at_extremity, nproc,nproc_x &
     &,nproc_y, neighbour, g_rows, g_row_length, g_datastart, mype      &
! field dimensions etc.
     &,row_length, rows, n_rows, land_field, model_levels, wet_levels   &
     &,bl_levels, st_levels, sm_levels, Ozone_levels, cloud_levels      &
     &,land_ice_points,soil_points,n_cca_lev,ntiles,salt_dim1,salt_dim2 &
     &,salt_dim3,tr_levels,tr_vars,co2_dim_len,co2_dim_row,co2_dim_lev  &
     &,n_arcl_species, n_arcl_compnts, i_arcl_compnts                   &
! model switches
     &,     model_domain, L_regular, L_SEC_VAR, L_EqT                   &
     &,     L_Rad_Step, L_Rad_Step_diag, L_Rad_Step_prog                &
     &,     L_Forcing, L_Timestep, L_Radiance, L_Wenyi                  &
     &,     LCAL360, L_MICROPHY, L_emcorr, L_climat_aerosol, Ltimer     &
     &,     L_gwd, L_use_ussp, l_taus_scale, l_fix_gwsatn, l_gwd_40km   &
     &,     l_ussp_opaque, sat_scheme, l_use_clearrh, L_ssice_albedo    &
     &,     L_RHCPT, L_Murk, l_murk_source, l_murk_bdry, L_MURK_RAD     &
     &,     L_DUST, L_SULPC_SO2, L_SULPC_NH3, L_SOOT, L_BIOMASS, L_OCFF &
! kdcorbin, 06/10 - added L_CO2_RADIATION flag
     &,l_co2_interactive,l_co2_radiation                                &
     &,Lflux_reset,L_clim_aero_hgt,L_HadGEM1_Clim_Aero                  &
     &,L_use_dust, L_use_sulphate_autoconv, L_auto_debias,              &
     &L_use_seasalt_autoconv,L_use_seasalt_indirect,L_use_seasalt_direct&
     &,     L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a          &
     &,     L_snow_albedo, L_ctile, L_radiation, L_rain                 &
     &,     L_INHOM_CLOUD, L_USE_BIOGENIC                               &
     &,     L_USE_SULPC_DIRECT, l_use_soot_direct, l_use_soot_indirect  &
     &,     l_use_soot_autoconv, l_use_bmass_direct                     &
     &,     l_use_bmass_indirect, l_use_bmass_autoconv                  &
     &,     l_use_ocff_direct, l_use_ocff_indirect, l_use_ocff_autoconv &
     &,     L_USE_SULPC_INDIRECT_SW, L_USE_SULPC_INDIRECT_LW            &
     &,     L_pc2, L_eacf, L_mr_physics1,l_cry_agg_dep,l_droplet_settle &
     &,     L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup                       &
     &,     L_USE_METHOX, L_rad_deg, L_TRIFFID                          &
     &,     l_use_stochem_ch4, L_it_melting, L_ukca, L_USE_ARCL         &
     &,     L_USE_SPEC_SEA, L_MOD_BARKER_ALBEDO, L_MOD_K_FLUX           &
     &,     L_SCVARY, L_VOLCTS                                          &
! model Parameters
     &,     RHcrit, cw_sea, cw_land                                     &
     &,     A_SW_segments, A_SW_seg_size, A_LW_segments, A_LW_seg_size  &
     &,     A_SW_radstep, A_LW_radstep, A_SW_radstep_diag               &
     &,     A_SW_radstep_prog, A_LW_radstep_diag, A_LW_radstep_prog     &
     &,     aero_bl_levels,INHOM_CLOUD_SW, INHOM_CLOUD_LW               &
     &,     DP_CORR_STRAT, DP_CORR_CONV, CO2_MMR, alpham,alphac,alphab  &
     &,     dtice,dt_bare,dalb_bare_wet,pen_rad_frac,SW_beta            &
     &,min_trop_level, max_trop_level, KAY_GWAVE, GWD_FRC, gwd_fsat     &
     &,O2MMR, N2OMMR                                                    &
     &,CH4MMR, C11MMR, C12MMR, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR&
! Pass down position of greenhouse gases in free_tracers
!                array, for chemical coupling
     &,     ngrgas ,grgas_addr                                          &
     &,     cloud_fraction_method,overlap_ice_liquid                    &
     &,     ice_fraction_method,ctt_weight,t_weight                     &
     &,     qsat_fixed,sub_cld,dbsdtbs_turb_0                           &
     &,     L_seq_mcr,L_autoc_3b,L_autolim_3b,L_autoconv_murk           &
     &,     ec_auto,N_drop_land                                         &
     &,     N_drop_sea,N_drop_land_cr,N_drop_sea_cr,Ntot_land, Ntot_sea &
     &,     x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic                    &
     &,     lsp_ei,lsp_fi,lsp_eic,lsp_fic                               &
! parameter for stochastic physics random parameters2
     &,     M_CI                                                        &
! Physical constants
     &,Lc,Lf,Cp, two_Omega, p_zero, kappa,R, g, Lapse, earth_radius, Pi &
! Vertical coordinate levels.
     &, r_rho_levels, r_theta_levels, r_at_u, r_at_v, eta_theta_levels  &
     &, eta_rho_levels, delta_lambda, delta_phi, lat_rot_NP, long_rot_NP&
! Time stepping information
     &,     timestep, radiation_timestep, radiation_tstep_diag          &
     &,     radiation_tstep_prog, I_year, I_day_number, I_hour          &
     &,     I_minute, I_second, timestep_number, PREVIOUS_TIME          &
!glr cable 1 line
     &, istep_cur                                                       &
! trig arrays
     &,     sin_theta_longitude, cos_theta_longitude                    &
     &,     FV_cos_theta_latitude, sin_theta_latitude                   &
! grid-dependent arrays
     &,     f3_at_u, true_longitude, true_latitude,                     &
! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &      STASHwork1,STASHwork2,STASHwork4,STASHwork6,STASHwork14     &
!
! Additional variables for SCM diagnostics
     &,     nSCMDpkgs, L_SCMDiags                                       &
!
! Data Fields.
     &,     THETA, Q, QCL, QCF, QCF2                                    &
     &,     QRAIN, QGRAUP, RHO, U, V                                    &
     &,     P, PSTAR, EXNER_RHO_LEVELS                                  &
     &,     EXNER_THETA_LEVELS, LAND                                    &
     &,     P_THETA_LEVELS, FRAC_LAND,frac_control                      &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &,land_index, RGRAIN_TILE,SNSOOT,NTML, CUMULUS                     &
     &,ICE_FRACTION, CCA, CCB, CCT, CCLWP, CCW_RAD, LCBASE              &
     &,TSTAR,TSTAR_LAND,TSTAR_SEA, TSTAR_SICE                           &
     &,SICE_ALB, LAND_ALB, SNODEP,SNODEP_SEA                            &
     &,ozone3D,SW_INCS                                                  &
     &, LW_INCS, DIRPAR                                                 &
     &,     O3_trop_level, O3_trop_height, T_trop_level, T_trop_height  &
     &,ZH, OROG_SD, OROG_GRAD_XX, OROG_GRAD_XY                          &
     &,OROG_GRAD_YY, CF_AREA, CF_BULK                                   &
     &,CF_LIQUID,CF_FROZEN,MURK_SOURCE,arcl                             &
     &,SOIL_ALB,LAI_PFT,SNODEP_TILE, FRAC_TYP                           &
     &,TSTAR_TILE,Z0_TILE,DOLR_FIELD,LW_DOWN,SW_TILE_RTS                &
     &,ES_SPACE_INTERP, RAD_MASK, CH4_STOCH                             &
     &, cos_zenith_angle,can_rad_mod                                    &
! IN/OUT
     &,theta_star, q_star, qcl_star, qcf_star, qcf2_star, qrain_star    &
     &,qgraup_star, cf_star, cfl_star, cff_star, R_u, R_v               &
     &,a_realhd(rh_energy_corr), NET_FLUX, NET_MFLUX                    &
     &,MURK, DUST_DIV1, DUST_DIV2                                       &
     &,DUST_DIV3, DUST_DIV4, DUST_DIV5                                  &
     &,DUST_DIV6, SO2, SO4_AITKEN                                       &
     &,SO4_ACCU, SO4_DISS, NH3, SOOT_NEW                                &
     &,SOOT_AGD, SOOT_CLD, BMASS_NEW                                    &
     &,BMASS_AGD,BMASS_CLD, OCFF_NEW, OCFF_AGD, OCFF_CLD, CO2, TRACER   &
     &,biogenic, A_INTHD(23)                                            &
!glr cable 
     &, surf_down_sw,alb_tile                                           &
     &, day,TILE_PTS,SM_LEVELS,TILE_INDEX                               &
     &, SNOW_TMP3L,SNOW_RHO1L,TSOIL_TILE,SNOW_FLG3L,LAND_ALBEDO_CABLE   &
     &, l_cable                                                         &
! OUT Fields
     &,     ls_rain, ls_snow, micro_tends, unscaled_dry_rho             &
     &,     photosynth_act_rad, rad_hr, surf_radflux, dolr, SW_tile     &
! Section information
     &,     maxsects,h_sect                                             &
! error information
     &,     ErrorStatus  )

! DEPENDS ON: timer
        If (Ltimer) Call timer('Atmos_Physics1',4)
      
      ! Deallocate the array of the aerosol climatology for NWP
      deallocate(arcl)
      
      ! --------------------------------------------------------------
      ! Convert back from  mixing ratios to spec. humidities if needed
      ! --------------------------------------------------------------
      If (l_mr_physics1) Then

        ! Allocate mix_star variables to hold mixing ratio increments
        ! from atmos_physics1

        allocate ( mix_v_star (1-offx:row_length+offx,                  &
     &                   1-offy:rows+offy, wet_levels) )
        allocate ( mix_cl_star(1-offx:row_length+offx,                  &
     &                   1-offy:rows+offy, wet_levels) )
        allocate ( mix_cf_star(1-offx:row_length+offx,                  &
     &                   1-offy:rows+offy, wet_levels) )

        If (L_mcr_qcf2) then
          allocate ( mix_cf2_star(1-offx:row_length+offx,               &
     &                   1-offy:rows+offy, wet_levels) )
        Else
          allocate ( mix_cf2_star(1,1,1) )
        End if

        If (L_mcr_qrain) then
          allocate ( mix_rain_star(1-offx:row_length+offx,              &
     &                   1-offy:rows+offy, wet_levels) )
        Else
          allocate ( mix_rain_star(1,1,1) )
        End if

        If (L_mcr_qgraup) then
          allocate ( mix_graup_star(1-offx:row_length+offx,             &
     &                   1-offy:rows+offy, wet_levels) )
        Else
          allocate ( mix_graup_star(1,1,1) )
        End if

        ! q_star currently holds mixing ratio increments from physics1
        ! but we wish it to hold specific humidity increments.
        ! Start by copying q_star into mix_v_star etc.

        Do k = 1, wet_levels
          Do j = 1-offy, rows+offy
            Do i = 1-offx, row_length+offx
              mix_v_star(i,j,k)  = q_star(i,j,k)     ! Vapour
              mix_cl_star(i,j,k) = qcl_star(i,j,k)   ! Liquid
              mix_cf_star(i,j,k) = qcf_star(i,j,k)   ! Ice
            End do
          End do
        End do

        If (L_mcr_qcf2) then  ! Ice2
          Do k = 1, wet_levels
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                mix_cf2_star(i,j,k) = qcf2_star(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qrain) then  ! Rain
          Do k = 1, wet_levels
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                mix_rain_star(i,j,k) = qrain_star(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qgraup) then  ! Graupel
          Do k = 1, wet_levels
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                mix_graup_star = qgraup_star(i,j,k)
              End do
            End do
          End do
        End if

        ! Now convert the mixing ratio increments (stored in
        ! mix_v_star etc.) to specific humidity increments
        ! and write these s.h. incs. back into the q_star variables.
        ! d1 variables currently contain the mixing ratios before
        ! atmos_physics1.
        ! q_store variables currently contain specific humidities
        ! before atmos_physics1.

! DEPENDS ON: calc_q_star
        call calc_q_star (row_length, rows, wet_levels,                 &
     &           halo_i, halo_j, offx, offy,                            &
     &           Q, QCL, QCF,                                           &
     &           QCF2, QRAIN, QGRAUP,                                   &
     &           mix_v_star, mix_cl_star, mix_cf_star,                  &
     &           mix_cf2_star, mix_rain_star, mix_graup_star,           &
     &           L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                 &
     &           q_store, qcl_store, qcf_store,                         &
     &           qcf2_store, qrain_store, qgraup_store,                 &
     &           q_star, qcl_star, qcf_star,                            &
     &           qcf2_star, qrain_star, qgraup_star                     &
     &                   )

        ! Copy contents of q_store (specific humidities before
        ! atmos_physics1) back into the d1(q) variables etc. Include
        ! the halo points.

        Do k = 1, wet_levels
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              q(i,j,k)   = q_store(i,j,k)    ! Vapour
              qcl(i,j,k) = q_store(i,j,k)    ! Liquid
              qcf(i,j,k) = q_store(i,j,k)    ! Ice
            End do
          End do
        End do

        If (L_mcr_qcf2) then  ! Ice2
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qcf2(i,j,k) = qcf2_store(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qrain) then  ! Rain
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qrain(i,j,k) = qrain_store(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qgraup) then  ! Graupel
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qgraup(i,j,k) = qgraup_store(i,j,k)
              End do
            End do
          End do
        End if

        ! At this point all the conversions are complete.
        ! d1(q) again contains the specific humidities before the
        ! atmos_physics1 call. q_star contains specific humidity
        ! increments from atmos_physics1. However, the actual
        ! calculations that were done within atmos_physics1 were done
        ! by using mixing ration variables.
        !
        !
        ! Finally, Deallocate the temporary variables

        deallocate(q_store)
        deallocate(qcl_store)
        deallocate(qcf_store)
        deallocate(qcf2_store)
        deallocate(qrain_store)
        deallocate(qgraup_store)
        deallocate(mix_v_star)
        deallocate(mix_cl_star)
        deallocate(mix_cf_star)
        deallocate(mix_cf2_star)
        deallocate(mix_rain_star)
        deallocate(mix_graup_star)

      ! -----------------------------------------------------------
      ! End of conversion to back to specific humidities
      ! -----------------------------------------------------------
      End if  ! l_mr_physics1

        ! Clear workspace for radiation
        DEALLOCATE ( ozone3D )
! Deallocate space for tropopause diagnostics
        DEALLOCATE (O3_trop_level)
        DEALLOCATE (O3_trop_height)
        DEALLOCATE (T_trop_level)
        DEALLOCATE (T_trop_height)

! Diagnostics STASHed for each section in Atmos_Physics1:
      If (L_radiation) then
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,1,STASHwork1,                              &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',4)
        DEALLOCATE (STASHwork1)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,2,STASHwork2,                              &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',4)
        DEALLOCATE (STASHwork2)
      Endif

      If (L_rain) then
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,4,STASHwork4,                              &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',4)
        DEALLOCATE (STASHwork4)
      Endif


      IF(L_gwd .or. l_use_ussp)then
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,6,STASHwork6,                              &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',4)
        DEALLOCATE (STASHwork6)
      endif

        IF(L_EMCORR) THEN
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,14,STASHwork14,                          &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork14)
        ENDIF   ! l_emcorr

      Else  ! L_physics =.false.
! initialise arrays that hold physics increments to zero, only needs
! doing at non-halo points
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              theta_star(i,j,k) = 0.0
            End Do
          End Do
        End Do
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              q_star(i,j,k) = 0.0
              qcl_star(i,j,k) = 0.0
              qcf_star(i,j,k) = 0.0
              cf_star(i,j,k) = 0.0
              cfl_star(i,j,k) = 0.0
              cff_star(i,j,k) = 0.0
            End Do
          End Do
        End Do
        ! Initialise additional microphysics variables if in use
        If (L_mcr_qcf2)   qcf2_star(1:row_length, 1:rows, :) = 0.0
        If (L_mcr_qrain)  qrain_star(1:row_length, 1:rows, :) = 0.0
        If (L_mcr_qgraup) qgraup_star(1:row_length, 1:rows, :) = 0.0
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u(i,j,k) = 0.0
            End Do
          End Do
        End Do
        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = 0.0
            End Do
          End Do
        End Do

      End If ! L_Physics

      If (L_tracer .and. ErrorStatus == 0) then

! store physics changes for use in Sl_tracer2

! DEPENDS ON: tr_set_phys
        call TR_Set_Phys(                                               &
                         super_array_size, super_tracer_phys1,          &
                         L_CO2_interactive, CO2,                        &
                         L_Murk_advect, murk,                           &
                         L_Soot, soot_new,                              &
                                 soot_agd,                              &
                                 soot_cld,                              &
                         L_SULPC_SO2, SO2,                              &
                                      SO4_aitken,                       &
                                      so4_accu,                         &
                                      so4_diss,                         &
                         L_sulpc_nh3, nh3,                              &
                         L_sulpc_dms, dms,                              &
                         L_dust, DUST_DIV1,                             &
                                 DUST_DIV2,                             &
                                 DUST_DIV3,                             &
                                 DUST_DIV4,                             &
                                 DUST_DIV5,                             &
                                 DUST_DIV6,                             &
                         L_biomass, bmass_new,                          &
                                    bmass_agd,                          &
                                    bmass_cld,                          &
                         L_ocff, ocff_new,                              &
                                 ocff_agd,                              &
                                 ocff_cld,                              &
                         L_USE_CARIOLLE, OZONE_TRACER,                  &
                         tracer_phys1, tracer, tracer_ukca,             &
                         row_length, rows,                              &
                         model_levels, tr_levels, tr_vars, tr_ukca,     &
                         offx, offy, model_domain,                      &
                         .true., halo_i, halo_j                         &
                                               )

      end if  ! L_tracer and ErrorStatus == 0
     
      IF (L_mix_ratio) Then
        allocate ( mix_v (1-halo_i:row_length+halo_i,                   &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        allocate ( mix_cl(1-halo_i:row_length+halo_i,                   &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        allocate ( mix_cf(1-halo_i:row_length+halo_i,                   &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        allocate ( mix_v_star (1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy,wet_levels) )
        allocate ( mix_cl_star (1-offx:row_length+offx,                 &
     &                         1-offy:rows+offy,wet_levels) )
        allocate ( mix_cf_star (1-offx:row_length+offx,                 &
     &                         1-offy:rows+offy,wet_levels) )
        IF(L_mcr_qcf2)then
          allocate ( mix_cf2      (1-halo_i:row_length+halo_i,          &
     &                             1-halo_j:rows+halo_j, wet_levels) )
          allocate ( mix_cf2_star (1-offx:row_length+offx,              &
     &                             1-offy:rows+offy,wet_levels) )
        else
          allocate ( mix_cf2      (1,1,1) )
          allocate ( mix_cf2_star (1,1,1) )
        endif
        IF(L_mcr_qrain)then
          allocate ( mix_rain     (1-halo_i:row_length+halo_i,          &
     &                             1-halo_j:rows+halo_j, wet_levels) )
          allocate ( mix_rain_star(1-offx:row_length+offx,              &
     &                             1-offy:rows+offy,wet_levels) )
        else
          allocate ( mix_rain      (1,1,1) )
          allocate ( mix_rain_star (1,1,1) )
        endif
        IF(L_mcr_qgraup)then
          allocate ( mix_graup     (1-halo_i:row_length+halo_i,         &
     &                             1-halo_j:rows+halo_j, wet_levels) )
          allocate ( mix_graup_star(1-offx:row_length+offx,             &
     &                             1-offy:rows+offy,wet_levels) )
        else
          allocate ( mix_graup      (1,1,1) )
          allocate ( mix_graup_star (1,1,1) )
        endif

!  convert q, qcl,qcf to mix_v, mix_cl,mix_cf
! DEPENDS ON: q_to_mix_halo
        call q_to_mix_halo (row_length, rows, wet_levels,               &
     &               halo_i, halo_j,                                    &
     &               Q, QCL, QCF,                                       &
     &               QCF2, QRAIN, QGRAUP,                               &
     &               L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,             &
     &               mix_v, mix_cl, mix_cf,                             &
     &               mix_cf2, mix_rain, mix_graup )

!  convert q, qcl,qcf _star to mix_v, mix_cl,mix_cf _star
!  q_star holds delta_q

! DEPENDS ON: calc_mix_star
        call calc_mix_star                                              &
     &                  (row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, offx, offy,                    &
     &                   Q, QCL, QCF,                                   &
     &                   QCF2, QRAIN, QGRAUP,                           &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   qcf2_star, qrain_star, qgraup_star,            &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   mix_v_star, mix_cl_star, mix_cf_star           &
     &                  ,mix_cf2_star, mix_rain_star, mix_graup_star    &
     &                   )

      Else

        allocate ( mix_v_star (1,1,1) )
        allocate ( mix_cl_star (1,1,1) )
        allocate ( mix_cf_star (1,1,1) )
        allocate ( mix_cf2_star (1,1,1) )
        allocate ( mix_rain_star (1,1,1) )
        allocate ( mix_graup_star (1,1,1) )

      End If    ! L_mix_ratio


! store physics changes for use in Sl_moist_conserve
      If (L_moist_nonhydro_conserve .and.                               &
     &      ErrorStatus == 0 )then
! At this point q_star holds the increment to q
        IF(l_mix_ratio)then
          allocate ( mix_v_phys1(1-offx:row_length+offx,                &
     &                     1-offy:rows+offy, wet_levels) )
          allocate ( mix_cl_phys1(1-offx:row_length+offx,               &
     &                       1-offy:rows+offy, wet_levels) )
          allocate ( mix_cf_phys1(1-offx:row_length+offx,               &
     &                       1-offy:rows+offy, wet_levels) )
          if(L_mcr_qcf2)then
            allocate ( mix_cf2_phys1(1-offx:row_length+offx,            &
     &                               1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_cf2_phys1(1,1,1) )
          endif
          if(L_mcr_qrain)then
            allocate ( mix_rain_phys1(1-offx:row_length+offx,           &
     &                                1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_rain_phys1(1,1,1) )
          endif
          if(L_mcr_qgraup)then
            allocate ( mix_graup_phys1(1-offx:row_length+offx,          &
     &                                 1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_graup_phys1(1,1,1) )
          endif
! DEPENDS ON: calc_mix_star
          call calc_mix_star                                            &
     &                  (row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, offx, offy,                    &
     &                   Q, QCL, QCF,                                   &
     &                   QCF2, QRAIN, QGRAUP,                           &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   qcf2_star, qrain_star, qgraup_star,            &
!    &                   .false.,.false.,.false.,
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   mix_v_phys1, mix_cl_phys1, mix_cf_phys1        &
     &                   ,mix_cf2_phys1, mix_rain_phys1, mix_graup_phys1&
     &                   )
          if(l_pc2)then
          allocate ( cf_phys1 (1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
          allocate ( cfl_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
          allocate ( cff_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                cf_phys1(i,j,k)  = cf_star(i,j,k)
                cfl_phys1(i,j,k) = cfl_star(i,j,k)
                cff_phys1(i,j,k) = cff_star(i,j,k)
              End Do
            End Do
          End Do
          else
            allocate ( cf_phys1 (1,1,1) )
            allocate ( cfl_phys1 (1,1,1) )
            allocate ( cff_phys1 (1,1,1) )
         endif
        else
          allocate ( q_phys1(1-offx:row_length+offx,                    &
     &                       1-offy:rows+offy, wet_levels) )
          allocate ( qcl_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
          allocate ( qcf_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
         if(l_pc2)then
          allocate ( cf_phys1(1-offx:row_length+offx,                   &
     &                        1-offy:rows+offy, wet_levels) )
          allocate ( cfl_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
          allocate ( cff_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
         else
           allocate ( cf_phys1 (1,1,1) )
           allocate ( cfl_phys1 (1,1,1) )
           allocate ( cff_phys1 (1,1,1) )
         endif
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_phys1(i,j,k) = q_star(i,j,k)
                qcl_phys1(i,j,k) = qcl_star(i,j,k)
                qcf_phys1(i,j,k) = qcf_star(i,j,k)
              End Do
            End Do
          End Do
          if(L_pc2)then
           Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  cf_phys1(i,j,k)  = cf_star(i,j,k)
                  cfl_phys1(i,j,k) = cfl_star(i,j,k)
                  cff_phys1(i,j,k) = cff_star(i,j,k)
                End Do
              End Do
            End Do
          endif
         if(L_mcr_qcf2)then
           allocate ( qcf2_phys1(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
           Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  qcf2_phys1(i,j,k)  = qcf2_star(i,j,k)
                End Do
              End Do
            End Do
         else
           allocate ( qcf2_phys1(1,1,1) )
         endif
         if(L_mcr_qrain)then
           allocate ( qrain_phys1(1-offx:row_length+offx,               &
     &                            1-offy:rows+offy, wet_levels) )
           Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  qrain_phys1(i,j,k)  = qrain_star(i,j,k)
                End Do
              End Do
            End Do
         else
           allocate ( qrain_phys1(1,1,1) )
         endif
         if(L_mcr_qgraup)then
           allocate ( qgraup_phys1(1-offx:row_length+offx,              &
     &                             1-offy:rows+offy, wet_levels) )
           Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  qgraup_phys1(i,j,k)  = qgraup_star(i,j,k)
                End Do
              End Do
            End Do
         else
           allocate ( qgraup_phys1(1,1,1) )
         endif
        endif       !L_mix_ratio
      Else If ( NumCycles > 1 .and. L_pc2                               &
     &                           .and. ErrorStatus == 0 ) Then
!
! When cycling and PC2 is used cf_star etc need to be reset (at the
! beginning of each new cycle) to the value they had when they
! exited Physics1(). The following arrays hold these values.
!
        allocate ( cf_phys1(1-offx:row_length+offx,                     &
     &                      1-offy:rows+offy, wet_levels) )
        allocate ( cfl_phys1(1-offx:row_length+offx,                    &
     &                         1-offy:rows+offy, wet_levels) )
        allocate ( cff_phys1(1-offx:row_length+offx,                    &
     &                         1-offy:rows+offy, wet_levels) )
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              cf_phys1(i,j,k)  = cf_star(i,j,k)
              cfl_phys1(i,j,k) = cfl_star(i,j,k)
              cff_phys1(i,j,k) = cff_star(i,j,k)
            End Do
          End Do
        End Do
      end if  !L_moist_nonhydro_conserve

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS1 Atmos_Phys1',6)
! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after Atmos_Phys1 ***'
        ENDIF ! L_print_pe .or. mype ==0
        R_w =0.0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      endif !  L_print_L2norms
! ---------------------------------------------------------------

!
! When cycling the following variables need to be reset (at the
! beginning of each new cycle) to the value they had when they
! exited Physics1().
!

      If ( NumCycles > 1 ) Then

        Allocate( R_u_phys1(row_length,rows,model_levels) )
        Allocate( R_v_phys1(row_length,n_rows,model_levels) )
        Allocate( thetastar_phys1(row_length,rows, model_levels) )
        Allocate( qstar_phys1(row_length,rows,wet_levels) )
        Allocate( qclstar_phys1(row_length,rows,wet_levels) )
        Allocate( qcfstar_phys1(row_length,rows,wet_levels) )
        If (L_mcr_qcf2)                                                 &
     &    Allocate( qcf2_star_phys1(row_length,rows,wet_levels) )
        If (L_mcr_qrain)                                                &
     &    Allocate( qrain_star_phys1(row_length,rows,wet_levels) )
        If (L_mcr_qgraup)                                               &
     &    Allocate( qgraup_star_phys1(row_length,rows,wet_levels) )
        If ( L_mix_ratio ) Then
          Allocate ( mix_v_star_phys1(row_length,rows,wet_levels) )
          Allocate ( mix_cl_star_phys1(row_length,rows,wet_levels) )
          Allocate ( mix_cf_star_phys1(row_length,rows,wet_levels) )
          If ( L_mcr_qcf2 )                                             &
     &      Allocate( mix_cf2_star_phys1(row_length,rows,wet_levels) )
          If ( L_mcr_qrain )                                            &
     &      Allocate( mix_rain_star_phys1(row_length,rows,wet_levels) )
          If ( L_mcr_qgraup )                                           &
     &      Allocate( mix_graup_star_phys1(row_length,rows,wet_levels) )
        End If
        Allocate( area_cld_frac_phys1(row_length,rows,wet_levels) )
        Allocate( bulk_cld_frac_phys1(row_length,rows,wet_levels) )
        Allocate( bulk_cld_liq_phys1(row_length,rows, wet_levels) )
        Allocate( bulk_cld_fr_phys1(row_length,rows, wet_levels) )
        Allocate( ti_phys1(row_length, rows, nice) )
        Allocate( zh_phys1(row_length, rows) )
        Allocate( z0msea_phys1(row_length, rows) )
        Allocate( cca_phys1 (row_length, rows, n_cca_lev) )
        Allocate( ccb_phys1 (row_length, rows ) )
        Allocate( cct_phys1 (row_length, rows ) )
        If ( L_ctile ) Then
          Allocate( T_LAND_CTILE_PHYS1(ROW_LENGTH,ROWS) )
          Allocate( T_SICE_CTILE_PHYS1(ROW_LENGTH,ROWS) )
        End If
        Allocate( t_surf_phys1(row_length, rows) )
        Allocate( t_sf_tile_phys1(land_field,ntiles) )
        Allocate( snow_tile_phys1(land_field,ntiles) )
        Allocate( dolr_phys1(row_length,rows) )

        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u_phys1(i,j,k)       = R_u(i,j,k)
              thetastar_phys1(i,j,k) = theta_star(i,j,k)
              qstar_phys1(i,j,k)     = q_star(i,j,k)
              qclstar_phys1(i,j,k)   = qcl_star(i,j,k)
              qcfstar_phys1(i,j,k)   = qcf_star(i,j,k)
            End Do
          End Do
        End Do

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v_phys1(i,j,k) = R_v(i,j,k)
            End Do
          End Do
        End Do

        Do k = wet_levels, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u_phys1(i,j,k)       = R_u(i,j,k)
              thetastar_phys1(i,j,k) = theta_star(i,j,k)
            End Do
          End Do
        Enddo

        If ( L_mcr_qcf2 ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qcf2_star_phys1(i,j,k) = qcf2_star(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mcr_qrain ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qrain_star_phys1(i,j,k) = qrain_star(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mcr_qgraup ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qgraup_star_phys1(i,j,k) = qgraup_star(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mix_ratio ) Then

          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                mix_v_star_phys1(i,j,k)  = mix_v_star(i,j,k)
                mix_cl_star_phys1(i,j,k) = mix_cl_star(i,j,k)
                mix_cf_star_phys1(i,j,k) = mix_cf_star(i,j,k)
              End do
            End do
          End do

          If ( L_mcr_qcf2 ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_cf2_star_phys1(i,j,k) = mix_cf2_star(i,j,k)
                End Do
              End Do
            End do
          End If

          If ( L_mcr_qrain ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_rain_star_phys1(i,j,k) = mix_rain_star(i,j,k)
                End Do
              End Do
            End do
          End If

          If ( L_mcr_qgraup ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_graup_star_phys1(i,j,k) = mix_graup_star(i,j,k)
                End Do
              End Do
            End do
          End If

        End If

        Do k = 1, wet_levels
          Do j = j_start, j_stop
!CDIR NODEP
            Do i = i_start, i_stop
              bulk_cld_frac_phys1(i,j,k) = cf_bulk(i,j,k)
              bulk_cld_liq_phys1(i,j,k)  = cf_liquid(i,j,k)
              bulk_cld_fr_phys1(i,j,k)   = cf_frozen(i,j,k)
              area_cld_frac_phys1(i,j,k) = cf_area(i,j,k)
            Enddo
          Enddo
        End do

        Do k=1, n_cca_lev
          Do j=j_start, j_stop
!CDIR NODEP
             Do i=i_start, i_stop
               cca_phys1(i,j,k) = cca(i,j,k)
             Enddo
          Enddo
        Enddo

        Do j=j_start, j_stop
!CDIR NODEP
          Do i=i_start, i_stop
            ccb_phys1(i,j)=ccb(i,j)
            cct_phys1(i,j)=cct(i,j)
          Enddo
        Enddo


! The following piece of code is intended to keep surface
! vars at the same timelevel at each cycle.
        If (nice == 1) Then
! Set sea ice catagory scheme D1 pointers as catogories are not pres
! if nice = 1.
          p_ti        => ti
        else
          p_ti        => ti_cat
        End If

        Do k=1, nice
          Do j = 1, rows
!CDIR NODEP
            Do i = 1, row_length
               ti_phys1(i,j,k) = p_ti(i,j,k)
            Enddo
          Enddo
        Enddo

        Do j = 1, rows
!CDIR NODEP
          Do i = 1, row_length
             z0msea_phys1(i,j) = z0(i,j)
             zh_phys1(i,j) = zh(i,j)
          Enddo
        Enddo

        Do j = 1, rows
!CDIR NODEP
          Do i = 1, row_length
            t_surf_phys1(i,j) = TSTAR(i,j)
            dolr_phys1(i,j)   = dolr(i,j)
          Enddo
        Enddo

        Do i=1, land_field
           GS1(i) = GS(i)
        Enddo

        If ( L_ctile ) Then
          Do j = 1, rows
!CDIR NODEP
            Do i = 1, row_length
              t_land_ctile_phys1(i,j) = TSTAR_LAND(i,j)
              t_sice_ctile_phys1(i,j) = TSTAR_SICE(i,j)
            Enddo
          Enddo
        End If

        Do j = 1, ntiles
!CDIR NODEP
          Do i = 1, land_field
            t_sf_tile_phys1(i,j) = TSTAR_TILE(i+(j-1)*land_field)
            snow_tile_phys1(i,j) = SNODEP_TILE(i+(j-1)*land_field)
          Enddo
        Enddo

      End If  ! Num_Cycles > 1

      If ( L_new_tdisc ) Then

        Allocate( theta_np1(1-offx:row_length+offx,                     &
     &               1-offy:rows+offy, model_levels))
        Allocate( rho_np1 (1-offx:row_length+offx,                      &
     &       1-offy:rows+offy, model_levels) )
        Allocate( u_np1(1-offx:row_length+offx,                         &
     &       1-offy:rows+offy, model_levels ) )
        Allocate( v_np1(1-offx:row_length+offx,                         &
     &       1-offy:n_rows+offy, model_levels) )
        Allocate( w_np1(1-offx:row_length+offx,                         &
     &                  1-offy:rows+offy, 0:model_levels) )

        If ( L_mix_ratio ) Then

          Allocate ( mix_v_np1(1-offx:row_length+offx,                  &
     &                          1-offy:rows+offy, wet_levels) )
          Allocate ( mix_cl_np1(1-offx:row_length+offx,                 &
     &                          1-offy:rows+offy, wet_levels) )
          Allocate ( mix_cf_np1(1-offx:row_length+offx,                 &
     &                          1-offy:rows+offy, wet_levels) )

          If ( L_mcr_qcf2 ) Then
            Allocate( mix_cf2_np1(1-offx:row_length+offx,               &
     &                          1-offy:rows+offy, wet_levels) )
          Else
            Allocate( mix_cf2_np1(1,1,1) )
          End If

          If ( L_mcr_qrain ) Then
            Allocate( mix_rain_np1(1-offx:row_length+offx,              &
     &                          1-offy:rows+offy, wet_levels) )
          Else
            Allocate( mix_rain_np1(1,1,1) )
          End If

          If ( L_mcr_qgraup ) Then
            Allocate( mix_graup_np1(1-offx:row_length+offx,             &
     &                          1-offy:rows+offy, wet_levels) )
          Else
            Allocate( mix_graup_np1(1,1,1) )
          End If

          Allocate( q_np1(1,1,1) )
          Allocate( qcl_np1(1,1,1) )
          Allocate( qcf_np1(1,1,1) )
          Allocate ( qcf2_np1(1,1,1) )
          Allocate ( qrain_np1(1,1,1) )
          Allocate ( qgraup_np1(1,1,1) )

        Else

          Allocate( q_np1(1-offx:row_length+offx,                   &
     &                    1-offy:rows+offy, wet_levels) )
          Allocate( qcl_np1(1-offx:row_length+offx,                 &
     &                      1-offy:rows+offy, wet_levels) )
          Allocate( qcf_np1(1-offx:row_length+offx,                 &
     &                      1-offy:rows+offy, wet_levels) )

          If ( L_mcr_qcf2 ) Then
            Allocate ( qcf2_np1(1-offx:row_length+offx,                 &
     &                          1-offy:rows+offy, wet_levels) )
          Else
            Allocate ( qcf2_np1(1,1,1) )
          End If

          If ( L_mcr_qrain ) Then
            Allocate ( qrain_np1(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
          Else
            Allocate ( qrain_np1(1,1,1) )
          End If

          If ( L_mcr_qgraup ) Then
            Allocate ( qgraup_np1(1-offx:row_length+offx,               &
     &                            1-offy:rows+offy, wet_levels) )
          Else
            Allocate ( qgraup_np1(1,1,1) )
          End If

          Allocate ( mix_v_np1(1,1,1) )
          Allocate ( mix_cl_np1(1,1,1) )
          Allocate ( mix_cf_np1(1,1,1) )
          Allocate( mix_cf2_np1(1,1,1) )
          Allocate( mix_rain_np1(1,1,1) )
          Allocate( mix_graup_np1(1,1,1) )

        End If

      Else

        Allocate( theta_np1(1,1,1) )
        Allocate( u_np1(1,1,1) )
        Allocate( v_np1(1,1,1) )
        Allocate( w_np1(1,1,1) )
        Allocate( rho_np1 (1,1,1) )
        Allocate( q_np1(1,1,1) )
        Allocate( qcl_np1(1,1,1) )
        Allocate( qcf_np1(1,1,1) )
        Allocate( qcf2_np1(1,1,1) )
        Allocate( qrain_np1(1,1,1) )
        Allocate( qgraup_np1(1,1,1) )
        Allocate( mix_v_np1(1,1,1) )
        Allocate( mix_cl_np1(1,1,1) )
        Allocate( mix_cf_np1(1,1,1) )
        Allocate( mix_cf2_np1(1,1,1) )
        Allocate( mix_rain_np1(1,1,1) )
        Allocate( mix_graup_np1(1,1,1) )

      End If  ! L_new_tdisc

! ----------------------------------------------------------------------
! Section 2. Calculate new time-level estimates to Theta, q, qcl, qcf
!            by calling semi-Lagrangian advection routine.
!            Estimated values are returned in the _star variables.
! ----------------------------------------------------------------------
        ALLOCATE (exner_prime( 1-offx:row_length+offx, 1-offy:rows+offy,&
     &                         model_levels) )

        exner_prime(:,:,:) = 0.0

!
! Iterative SISL: Cycle through dynamics-physics to enable trajectory
! calculations from interpolated winds and utilize improved time
! discretization if requested.
!
      Do CycleNo = 1, NumCycles

! Restore phys1 variables to be used as predictors.
      If ( CycleNo > 1 ) Then

! reset weights after the first cycle
        alpha1 = alpha_1_2
        alpha2 = alpha_2_2
        alpha3 = alpha_3_2
        alpha4 = alpha_4_2

        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u(i,j,k)        = R_u_phys1(i,j,k)
              theta_star(i,j,k) = thetastar_phys1(i,j,k)
              q_star(i,j,k)     = qstar_phys1(i,j,k)
              qcl_star(i,j,k)   = qclstar_phys1(i,j,k)
              qcf_star(i,j,k)   = qcfstar_phys1(i,j,k)
            End Do
          End Do
        Enddo

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = R_v_phys1(i,j,k)
            End Do
          End Do
        End Do

        Do k = wet_levels, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u(i,j,k)        = R_u_phys1(i,j,k)
              theta_star(i,j,k) = thetastar_phys1(i,j,k)
            End Do
          End Do
        End Do

        If ( L_mcr_qcf2 ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qcf2_star(i,j,k) = qcf2_star_phys1(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mcr_qrain ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qrain_star(i,j,k) = qrain_star_phys1(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mcr_qgraup ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qgraup_star(i,j,k) = qgraup_star_phys1(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mix_ratio ) Then

          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                mix_v_star(i,j,k)  = mix_v_star_phys1(i,j,k)
                mix_cl_star(i,j,k) = mix_cl_star_phys1(i,j,k)
                mix_cf_star(i,j,k) = mix_cf_star_phys1(i,j,k)
              End do
            End do
          End do

          If ( L_mcr_qcf2 ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_cf2_star(i,j,k) = mix_cf2_star_phys1(i,j,k)
                End Do
              End Do
            End do
          End If

          If ( L_mcr_qrain ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_rain_star(i,j,k) = mix_rain_star_phys1(i,j,k)
                End Do
              End Do
            End do
          End If

          If ( L_mcr_qgraup ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_graup_star(i,j,k) = mix_graup_star_phys1(i,j,k)
                End Do
              End Do
            End do
          End If

        End If ! L_mix_ratio

        Do k = 1, wet_levels
          Do j = j_start, j_stop
!CDIR NODEP
            Do i = i_start, i_stop
              cf_bulk(i,j,k)   = bulk_cld_frac_phys1(i,j,k)
              cf_liquid(i,j,k) = bulk_cld_liq_phys1(i,j,k)
              cf_frozen(i,j,k) = bulk_cld_fr_phys1(i,j,k)
              cf_area(i,j,k)   = area_cld_frac_phys1(i,j,k)
            Enddo
          Enddo
        End do

        If ( L_pc2 ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                cf_star(i,j,k)  = cf_phys1(i,j,k)
                cfl_star(i,j,k) = cfl_phys1(i,j,k)
                cff_star(i,j,k) = cff_phys1(i,j,k)
              End Do
            End Do
          End Do
        End If


        Do k=1, n_cca_lev
          Do j=j_start, j_stop
!CDIR NODEP
             Do i=i_start, i_stop
               cca(i,j,k) = cca_phys1(i,j,k)
             Enddo
          Enddo
        Enddo

        Do j=j_start, j_stop
!CDIR NODEP
          Do i=i_start, i_stop
            ccb(i,j) = ccb_phys1(i,j)
            cct(i,j) = cct_phys1(i,j)
          Enddo
        Enddo

        Do k=1, nice
          Do j = 1, rows
!CDIR NODEP
            Do i = 1, row_length
               p_ti(i,j,k) = ti_phys1(i,j,k)
            Enddo
          Enddo
        Enddo

        Do j = 1, rows
!CDIR NODEP
          Do i = 1, row_length
            z0(i,j)  = z0msea_phys1(i,j)
            zh(i,j)  = zh_phys1(i,j)
          Enddo
        Enddo

        Do i=1, land_field
           GS(i) = GS1(i)
        Enddo

        Do j = 1, rows
!CDIR NODEP
          Do i = 1, row_length
            TSTAR(i,j) = t_surf_phys1(i,j)
            dolr(i,j) = dolr_phys1(i,j)
          Enddo
        Enddo

        If ( L_ctile ) Then
          Do j = 1, rows
!CDIR NODEP
            Do i = 1, row_length
              TSTAR_LAND(i,j) = t_land_ctile_phys1(i,j)
              TSTAR_SICE(i,j) = t_sice_ctile_phys1(i,j)
            Enddo
          Enddo
        End If

        Do j = 1, ntiles
!CDIR NODEP
          Do i = 1, land_field
            TSTAR_TILE(i+(j-1)*land_field)  = t_sf_tile_phys1(i,j)
            SNODEP_TILE(i+(j-1)*land_field) = snow_tile_phys1(i,j)
          Enddo
        Enddo

      Else

! umui user defined alphas at first cycle
        alpha1 = alpha_1
        alpha2 = alpha_2
        alpha3 = alpha_3
        alpha4 = alpha_4

      End If

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS2 S-L Advection',5)

! Note: The departure points which are passed into sl_full_wind are
! returned in depart_lambda for lambda, depart_phi for phi and
! depart_r_w for r. depart_r_theta is passed into sl_tracer1.

      If (ErrorStatus  ==  0 ) Then

! Apply diagnostics at last cycle only.
       If ( CycleNo == NumCycles ) Then

! Save current values to form diagnostics increments over advection
        IF( sf(185,12) ) THEN
          ALLOCATE ( u_incr_diagnostic(row_length,rows,model_levels) )
! Hold u increment (no halos needed)
          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                u_incr_diagnostic(i,j,k) =R_u(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                   ! on STASHflag

        IF( sf(186,12) ) THEN
          ALLOCATE ( v_incr_diagnostic(row_length,n_rows,model_levels) )
! Hold v increment (no halos needed)
          DO k=1,model_levels
            DO j=1,n_rows
              DO i=1,row_length
                v_incr_diagnostic(i,j,k) =R_v(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(181,12) ) THEN
          ALLOCATE ( T_incr_diagnostic(row_length,rows,model_levels) )
! note at this point theta_star holds all increments to theta
          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                T_incr_diagnostic(i,j,k) =theta_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(182,12) ) THEN
          ALLOCATE ( q_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point q_star holds all increments to q
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                q_incr_diagnostic(i,j,k) =q_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(183,12) ) THEN
          ALLOCATE ( qcl_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point qcl_star holds all increments to qcl
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                qcl_incr_diagnostic(i,j,k) =qcl_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(184,12) ) THEN
          ALLOCATE ( qcf_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point qcf_star holds all increments to qcf
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                qcf_incr_diagnostic(i,j,k) =qcf_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(192,12) ) THEN
          ALLOCATE ( cf_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point cf_star holds all increments to cf
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                cf_incr_diagnostic(i,j,k) =cf_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(193,12) ) THEN
          ALLOCATE ( cfl_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point cfl_star holds all increments to cfl
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                cfl_incr_diagnostic(i,j,k) =cfl_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(194,12) ) THEN
          ALLOCATE ( cff_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point cff_star holds all increments to cff
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                cff_incr_diagnostic(i,j,k) =cff_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

       End If ! CycleNo == NumCycles

!      IF( sf(187,12) ) THEN    ! w increments
! Note no need to hold R_w before semi-lagrangian advection as unset.

        If ( CycleNo == 1 ) Then
        ALLOCATE ( depart_lambda(row_length,rows,model_levels) )
        ALLOCATE ( depart_phi(row_length,rows,model_levels) )
        ALLOCATE ( depart_r_w(row_length,rows,model_levels) )
        ALLOCATE ( depart_r_theta(row_length,rows,model_levels) )
        End If

! DEPENDS ON: timer
        If (Ltimer) Call timer ('SL_Thermo',3)

! DEPENDS ON: ni_sl_thermo
        Call NI_SL_Thermo(                                              &
     &                    moisture_array_size,                          &
     &                    theta, q, qcl, qcf, qcf2, qrain, qgraup,      &
     &                    mix_v, mix_cl, mix_cf,                        &
     &                    mix_cf2, mix_rain, mix_graup,                 &
     &                    cf_bulk, cf_liquid, cf_frozen,                &
     &                    q_star, qcl_star, qcf_star,                   &
     &                    qcf2_star, qrain_star, qgraup_star,           &
     &                    mix_v_star, mix_cf_star, mix_cl_star,         &
     &                    mix_cf2_star, mix_rain_star, mix_graup_star,  &
     &                    cf_star, cfl_star, cff_star,                  &
     &                    exner_star, theta_star, theta_np1,            &
     &                    w, w_adv, u_adv, v_adv,                       &
     &                    exner_theta_levels,                           &
     &                    pstar, p, p_theta_levels, rho,                &
     &             eta_rho_levels, eta_theta_levels, r_rho_levels,      &
     &             r_theta_levels, row_length, rows, n_rows,            &
     &                    model_levels, wet_levels,                     &
     &                    alpha_2, check_bottom_levels,                 &
     &                    interp_vertical_search_tol,                   &
     &                    first_constant_r_rho_level,                   &
     &                    delta_lambda, delta_phi,                      &
     &             glambda_p, phi_p, glambda_u, phi_v,                  &
     &             gdlambda_p, dphi_p, gdlambda_u, dphi_v, grecip_dlamp,&
     &             recip_dphip, grecip_dlamu, recip_dphiv,              &
     &             wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v,        &
     &             lambda_p_rm, lambda_p_rp, lambda_u_rm, lambda_u_rp,  &
     &             phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,              &
     &           recip_lambda_p_m, recip_lambda_p_0, recip_lambda_p_p,  &
     &           recip_lambda_p_p2, recip_lambda_u_m, recip_lambda_u_0, &
     &           recip_lambda_u_p, recip_lambda_u_p2, recip_phi_p_m,    &
     &           recip_phi_p_0, recip_phi_p_p, recip_phi_p_p2,          &
     &           recip_phi_v_m, recip_phi_v_0, recip_phi_v_p,           &
     &           recip_phi_v_p2, Base_lambda, base_phi, lambda_p_end,   &
     &           phi_p_end, dlambda_p_end, dphi_p_end, dphi_v_end,      &
     &           recip_dlam, recip_dphi, max_look,                      &
     &           look_lam, look_phi, halo_lam, halo_phi,                &
     &                 timestep, FV_cos_theta_latitude,                 &
     &                 cos_theta_latitude, sec_theta_latitude,          &
     &                 sin_theta_latitude, tan_theta_latitude,          &
     &                 cos_v_latitude, sec_v_latitude,                  &
     &                 sin_v_latitude, tan_v_latitude,                  &
     &                 sin_theta_longitude, cos_theta_longitude,        &
     &                 LAM_max_cfl, THETA_LBC, n_rims_to_do,            &
     &                 RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,           &
     &                 LENRIMA(fld_type_p,halo_type_extended,           &
     &                         rima_type_norm),                         &
     &                 LBC_SIZEA(1,fld_type_p,halo_type_extended,       &
     &                           rima_type_norm),                       &
     &                 LBC_STARTA(1,fld_type_p,halo_type_extended,      &
     &                            rima_type_norm),                      &
     &                 r_at_u, r_at_v,                                  &
     &                 mype, nproc, nproc_x, nproc_y,                   &
     &                 offx, offy, halo_i, halo_j, datastart,           &
     &                 g_i_pe, at_extremity,                            &
     &                 global_row_length, global_rows,                  &
     &                 gc_proc_row_group, gc_proc_col_group,            &
     &                 Pi, Depart_scheme, Depart_order,                 &
     &                 high_order_scheme(Theta_SL),                     &
     &                 monotone_scheme(Theta_SL),                       &
     &                 high_order_scheme(moist_SL),                     &
     &                 monotone_scheme(moist_SL),                       &
     &                 L_Ritchie_high, L_Ritchie_mono,                  &
     &                 Ritchie_high_order_scheme,                       &
     &                 Ritchie_monotone_scheme,                         &
     &                 model_domain, L_high(Theta_SL),                  &
     &                 L_mono(Theta_SL), thmono_levels,                 &
     &                 L_high(moist_SL), L_mono(moist_SL),              &
     &                 L_conserv(moist_SL), L_pc2,                      &
     &                 L_2d_sl_geometry, L_mix_ratio,                   &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 L_regular, L_sl_halo_reprod,                     &
     &                 L_fint_theta, L_lbc_new,                         &
     &                 L_new_tdisc, CycleNo,                            &
     &                 depart_r_theta,                                  &
     &                 depart_lambda, depart_phi, depart_r_w,           &
     &                 ErrorStatus )

        If (.not. L_moist_nonhydro_conserve .and.                       &
     &      .not. L_tracer .and. CycleNo == NumCycles ) Then
          DEALLOCATE (depart_r_theta)
        endif

! DEPENDS ON: timer
        If (Ltimer) Call timer ('SL_Thermo',4)

      End If        !     ErrorStatus  ==  0

! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after sl_thermo  ***'
        ENDIF ! L_print_pe .or. mype ==0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      endif !  L_print_L2norms

! ----------------------------------------------------------------------
! Section 2.1 Calculate advective Momentum increments.
! ----------------------------------------------------------------------
      If (ErrorStatus  ==  0 ) Then

! DEPENDS ON: timer
        If (Ltimer) Call timer ('SL_Full_wind',3)

! DEPENDS ON: ni_sl_full_wind
        Call NI_SL_Full_wind(                                           &
     &                      u, u_np1, v, v_np1,                         &
     &                      w, w_np1,                                   &
     &                      u_adv, v_adv, w_adv,                        &
     &                      theta, theta_np1,                           &
     &                      exner_rho_levels,                           &
     &                      q, qcl, qcf,                                &
     &                      qcf2, qrain, qgraup,                        &
     &                      q_np1, qcl_np1, qcf_np1, qcf2_np1,          &
     &                      qrain_np1, qgraup_np1,                      &
     &                      mix_v, mix_cl, mix_cf,                      &
     &                      mix_v_np1, mix_cl_np1, mix_cf_np1,          &
     &                      mix_cf2, mix_rain, mix_graup,               &
     &                      mix_cf2_np1, mix_rain_np1, mix_graup_np1,   &
     &                      L_mix_ratio,                                &
     &                      L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
     &                      depart_lambda, depart_phi, depart_r_w,      &
     &                      r_theta_levels, r_rho_levels,               &
     &                      eta_theta_levels, eta_rho_levels,           &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      sin_theta_latitude, cos_v_latitude,         &
     &                      sec_v_latitude, sin_v_latitude,             &
     &                      tan_theta_latitude, tan_v_latitude,         &
     &                      cos_theta_longitude,                        &
     &                      sin_theta_longitude,                        &
     &                      f1_at_v, f2_at_u, f3_at_u, f3_at_v,         &
     &                      delta_lambda, delta_phi, timestep,          &
     &                      glambda_p, phi_p, glambda_u, phi_v,         &
     &                      gdlambda_p, dphi_p, gdlambda_u, dphi_v,     &
     &            grecip_dlamp, recip_dphip, grecip_dlamu, recip_dphiv, &
     &            wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v,         &
     &            lambda_p_rm, lambda_p_rp, lambda_u_rm, lambda_u_rp,   &
     &            phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,               &
     &            recip_lambda_p_m, recip_lambda_p_0, recip_lambda_p_p, &
     &            recip_lambda_p_p2, recip_lambda_u_m, recip_lambda_u_0,&
     &            recip_lambda_u_p, recip_lambda_u_p2, recip_phi_p_m,   &
     &            recip_phi_p_0, recip_phi_p_p, recip_phi_p_p2,         &
     &            recip_phi_v_m, recip_phi_v_0, recip_phi_v_p,          &
     &            recip_phi_v_p2, base_lambda, base_phi, lambda_p_end,  &
     &            phi_p_end, dlambda_p_end, dphi_p_end, dphi_v_end,     &
     &            recip_dlam, recip_dphi, max_look,                     &
     &            look_lam, look_phi, halo_lam, halo_phi,               &
     &                      alpha_3, alpha_4, LAM_max_cfl,              &
     &                      n_Y_arrays, n_Yw_arrays,                    &
     &                      n_Yd_arrays, n_Ydw_arrays,                  &
     &                      U_LBC,V_LBC,W_LBC,                          &
     &                      LENRIMA(1,1,rima_type_norm),                &
     &                      LBC_SIZEA(1,1,1,rima_type_norm),            &
     &                      LBC_STARTA(1,1,1,rima_type_norm),           &
     &                      RIMWIDTHA(rima_type_norm), n_rims_to_do,    &
     &                      Pi, Cp, epsilon, g,                         &
     &                      model_domain, row_length, rows, n_rows,     &
     &                      model_levels, wet_levels,                   &
     &                      Depart_scheme, Depart_order,                &
     &                      high_order_scheme(Wind_SL),                 &
     &                      monotone_scheme(Wind_SL),                   &
     &                      L_trivial_trigs,                            &
     &                      L_high(Wind_SL), L_mono(Wind_SL),           &
     &                      L_conserv(Wind_SL),                         &
     &                      L_Ritchie_high,                             &
     &                      L_Ritchie_mono,                             &
     &                      Ritchie_high_order_scheme,                  &
     &                      Ritchie_monotone_scheme,                    &
     &                      first_constant_r_rho_level,                 &
     &                      check_bottom_levels,                        &
     &                      interp_vertical_search_tol,                 &
     &                      r_at_u, r_at_v,                             &
     &                      mype, nproc, nproc_x, nproc_y,              &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      global_row_length, global_rows,             &
     &                      datastart, at_extremity, g_i_pe,            &
     &                      gc_proc_row_group, gc_proc_col_group,       &
     &                      L_2d_sl_geometry, L_sl_halo_reprod,         &
     &                      L_free_slip, L_regular,                     &
     &                      L_qwaterload, L_interp_depart,              &
     &                      L_new_tdisc, CycleNo,                       &
     &                      R_u, R_v, R_w, ErrorStatus )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('SL_Full_wind',4)

! ---------------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after SL_Full_wind  ***'
        ENDIF ! L_print_pe .or. mype ==0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      endif !  L_print_L2norms
! ---------------------------------------------------------------

        IF ( .Not.SF(0,12) .and. CycleNo == NumCycles ) Then
          DEALLOCATE (depart_r_w)
          If (.not. L_moist_nonhydro_conserve .and.                     &
     &        .not. L_tracer )then
            DEALLOCATE (depart_phi)
            DEALLOCATE (depart_lambda)
          endif
        endif
      End If         !  ErrorStatus  ==  0

      IF(L_mix_ratio)then
!  convert mix_v, mix_cl,mix_cf to q, qcl,qcf

! DEPENDS ON: mix_to_q
        call mix_to_q                                                   &
     &                  (row_length, rows, wet_levels,                  &
     &                   offx, offy,                                    &
     &                   mix_v_star, mix_cl_star, mix_cf_star,          &
     &                   mix_cf2_star, mix_rain_star, mix_graup_star,   &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   qcf2_star, qrain_star, qgraup_star             &
     &                   )

      endif         !     L_mix_ratio

      If ( L_tracer .and. CycleNo == NumCycles ) Then

! DEPENDS ON: timer
        If (Ltimer) CALL TIMER('SL_tracer1',3)

! DEPENDS ON: sl_tracer1
        Call SL_tracer1(                                                &
     &                 super_array_size, eta_theta_levels,              &
     &                 r_rho_levels, r_theta_levels,                    &
     &                 PSTAR, P, P_THETA_LEVELS,                        &
     &                 RHO, l_tracer1_non_hydro,                        &
     &                 row_length, rows, n_rows, model_levels,          &
     &                 delta_lambda, delta_phi,                         &
     &                 glambda_p, phi_p, gdlambda_u, dphi_v,            &
     &                 grecip_dlamp, recip_dphip,                       &
     &                 lambda_p_rm, lambda_p_rp, phi_p_rm, phi_p_rp,    &
     &                 recip_lambda_p_m, recip_lambda_p_0,              &
     &                 recip_lambda_p_p, recip_lambda_p_p2,             &
     &                 recip_phi_p_m, recip_phi_p_0,                    &
     &                 recip_phi_p_p, recip_phi_p_p2,                   &
     &                 base_lambda, base_phi,                           &
     &                 recip_dlam, recip_dphi, max_look,                &
     &                 look_lam, look_phi, halo_lam, halo_phi,          &
     &                 FV_cos_theta_latitude,                           &
     &                 mype, nproc, nproc_x, nproc_y,                   &
     &                 halo_i, halo_j, datastart,                       &
     &                 g_i_pe, at_extremity,                            &
     &                 global_row_length, global_rows,                  &
     &                 gc_proc_row_group,                               &
     &                 gc_proc_col_group, offx, offy,                   &
     &                 L_regular, L_sl_halo_reprod,                     &
     &                 high_order_scheme(moist_SL),                     &
     &                 monotone_scheme(moist_SL),                       &
     &                 model_domain, L_high(moist_SL),                  &
     &                 L_mono(moist_SL),                                &
     &                 L_conserve_tracers,                              &
     &                 check_bottom_levels,                             &
     &                 interp_vertical_search_tol,                      &
     &                 first_constant_r_rho_level,                      &
     &                 depart_lambda, depart_phi, depart_r_theta,       &

!  d1(jxxx) holds time level n value plus physics1 increment
     &                 CO2, L_CO2_interactive,                          &
     &                 MURK, L_Murk_advect,                             &
     &                 DUST_DIV1,DUST_DIV2,                             &
     &                 DUST_DIV3,DUST_DIV4,                             &
     &                 DUST_DIV5,DUST_DIV6, L_DUST,                     &
     &                 SOOT_NEW, SOOT_AGD,                              &
     &                 SOOT_CLD, L_soot,                                &
     &                 BMASS_NEW, BMASS_AGD,                            &
     &                 BMASS_CLD, L_biomass,                            &
     &                 OCFF_NEW, OCFF_AGD, OCFF_CLD, L_OCFF,            &
     &                 SO2, SO4_AITKEN,                                 &
     &                 SO4_ACCU,                                        &
     &                 SO4_DISS, NH3, DMS,                              &
     &                 L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,           &
     &                 TRACER, tr_levels, tr_vars,                      &
     &                 tracer_ukca, tr_ukca,                            &
     &                 L_USE_CARIOLLE,    OZONE_TRACER,                 &
     &                 ErrorStatus)

! DEPENDS ON: timer
        If (Ltimer) CALL TIMER('SL_tracer1',4)

! ----------------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after SL_tracer1  ***'
        ENDIF ! L_print_pe .or. mype ==0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      endif !  L_print_L2norms
! ---------------------------------------------------------------

      end if  ! L_tracer

! ----------------------------------------------------------------------
! Section 2.21 Calculation of coefficients in turbulence scheme
! ----------------------------------------------------------------------
      If (L_subfilter_horiz .or. L_subfilter_vert) Then

        If (bl_levels  /=  model_levels - 1) then
          write(6,*)'BL_LEVELS =',bl_levels                             &
     &,             'MODEL_LEVELS =',model_levels
          ErrorStatus=123
! DEPENDS ON: ereport
          Call Ereport("ATM_STEP", ErrorStatus,                         &
     &        "The number of boundary layer levels must equal " //      &
     &        "the number of model levels minus one when the " //       &
     &        "3D subgrid turbulence scheme is selected.")
        End If

        If (diff_factor  <   0.0 .or. diff_factor  >   1.0) then
          write(6,*)'DIFF_FACTOR =',diff_factor
          ErrorStatus=123
! DEPENDS ON: ereport
          Call Ereport("ATM_STEP", ErrorStatus,                         &
     &        "DIFF_FACTOR must have a value greater than " //          &
     &        "zero and less than or equal to 1.0 so that the "//       &
     &        "numerical stability is maintained")
        End If

        If (mix_factor  <=  0.0 .or. mix_factor  >   1.0) then
          write(6,*)'MIX_FACTOR =',mix_factor
          ErrorStatus=123
! DEPENDS ON: ereport
          Call Ereport("ATM_STEP", ErrorStatus,                         &
     &        "MIX_FACTOR should have a value greater or equal " //     &
     &        "to zero and less than or equal to 1.0 ")
        End If

        If (turb_startlev_horiz > turb_endlev_horiz .OR.                &
     &      turb_startlev_vert > turb_endlev_vert) Then
          ErrorStatus=123
! DEPENDS ON: ereport
          Call Ereport("ATM_STEP", ErrorStatus,                         &
     &        "The start level for the turbulence scheme is " //        &
     &        "greater than the end level! ")
        End If

! The levels over which the turbulence scheme acts must be
! between 2 and model_levels-1

        If (turb_startlev_vert < 2) turb_startlev_vert = 2
        If (turb_endlev_vert > model_levels - 1)                        &
     &                           turb_endlev_vert = model_levels - 1

        If (turb_startlev_horiz < 2) turb_startlev_horiz = 2
        If (turb_endlev_horiz > model_levels - 1)                       &
     &                           turb_endlev_horiz = model_levels - 1

        If ( CycleNo == 1 ) Then
          ALLOCATE (visc_m(1-halo_i:row_length+halo_i                   &
     &                           , 1-halo_j:rows+halo_j, model_levels) )
          ALLOCATE (visc_h(1-halo_i:row_length+halo_i                   &
     &                           , 1-halo_j:rows+halo_j, model_levels) )
          ALLOCATE (RNEUTML(1:row_length, 1:rows,model_levels))
          ALLOCATE (shear(1-offx:row_length+offx, 1-offy:rows+offy      &
     &                                            , model_levels))
        End If

        Do k = 1, model_levels
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              visc_m(i,j,k)= 0.0
              visc_h(i,j,k)= 0.0
            End Do
          End Do
        End Do

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              shear(i,j,k)= 0.0
              RNEUTML(i,j,k) = 0.0
            End Do
          End Do
        End Do

        If (L_subfilter_vert) Then

          ALLOCATE (visc_BL_m(1:row_length, 1:rows, bl_levels))

          Do k = 1, bl_levels
            Do j = 1, rows
              Do i = 1, row_length
                visc_BL_m(i,j,k) = 0.0
              End Do
            End Do
          End Do

        Else

          ALLOCATE (visc_BL_m(1,1,1) )

        End If    ! L_subfilter_vert

!
! Calculate  lambda^2*S in TURB_Smagorinsky
!
! DEPENDS ON: turb_smagorinsky
        CALL TURB_Smagorinsky(                                          &
     &                     rows, row_length, n_rows                     &
     &,                    model_levels                                 &
     &,                    r_theta_levels, r_rho_levels                 &
     &,                    U, V, W                                      &
     &,                    visc_m, shear                                &
     &,                    Z0, RNEUTML, timestep                        &
     &,                    diff_factor, mix_factor, max_diff            &
     &,                    cos_theta_latitude                           &
     &,                    cos_v_latitude                               &
     &,                    delta_lambda, delta_phi)


        If (L_subfilter_vert) Then

          Do k = 2, bl_levels
            Do j = 1, rows
              Do i = 1, row_length
                visc_BL_m(i,j,k)= visc_m(i,j,k-1)
              End Do
            End Do
          End Do

        End If ! L_subfilter_vert

      Else

        If ( CycleNo == 1 ) Then
          ALLOCATE (visc_m(1,1,1) )
          ALLOCATE (visc_h(1,1,1) )
          ALLOCATE (visc_BL_m(1,1,1) )
          ALLOCATE (RNEUTML(1,1,1) )
          ALLOCATE (shear(1,1,1) )
        End If

      End If  !L_subfilter_horiz or L_subfilter_vert

! ---------------------------------------------------------------
! Section 2.3 Diagnostics at end of advection
! ----------------------------------------------------------------------
! Apply diagnostics at final cycle only
      If ( CycleNo == NumCycles ) Then

! DEPENDS ON: timer
      If (Ltimer) Call timer ('Diag_adv',3)

! section 12: 'dynamics advection' based quantities
      IF(      SF(0,12)                                                 &
                                ! Diagnostics required for this section
     &   .AND. ErrorStatus == 0) THEN

! Allocate diagnostic space for STASH
        ALLOCATE (STASHwork12(STASH_maxlen(12,A_im)))

! DEPENDS ON: diagnostics_adv
        CALL Diagnostics_adv(                                           &
     &            row_length, rows, n_rows, model_levels, wet_levels,   &
! primary wind fields:
     &            U, V,                                                 &
     &            THETA, Q, QCL, QCF,                                   &
     &            CF_BULK, CF_LIQUID, CF_FROZEN,                        &
! wind field increments after advection:
     &            R_u, R_v, R_w,                                        &
! wind field increments before advection (on stashflag):
     &            u_incr_diagnostic, v_incr_diagnostic,                 &
     &            T_incr_diagnostic, q_incr_diagnostic,                 &
     &            qcl_incr_diagnostic, qcf_incr_diagnostic,             &
     &            cf_incr_diagnostic, cfl_incr_diagnostic,              &
     &            cff_incr_diagnostic,                                  &
     &            theta_star, q_star, qcl_star, qcf_star,               &
     &            cf_star, cfl_star, cff_star,                          &
     &            EXNER_THETA_LEVELS,                                   &
! Departure points for w
     &            depart_lambda, depart_phi, depart_r_w,                &
     &            r_theta_levels,                                       &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &            STASHwork12)

        DEALLOCATE (depart_r_w)
        If (.not. L_moist_nonhydro_conserve .and.                       &
     &      .not. L_tracer )then
          DEALLOCATE (depart_phi)
          DEALLOCATE (depart_lambda)
        endif

! Tidy allocatable arrays
        IF( sf(185,12) ) THEN
          DEALLOCATE ( u_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(186,12) ) THEN
          DEALLOCATE ( v_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(181,12) ) THEN
          DEALLOCATE ( T_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(182,12) ) THEN
          DEALLOCATE ( q_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(183,12) ) THEN
          DEALLOCATE ( qcl_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(184,12) ) THEN
          DEALLOCATE ( qcf_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(192,12) ) THEN
          DEALLOCATE ( cf_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(193,12) ) THEN
          DEALLOCATE ( cfl_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(194,12) ) THEN
          DEALLOCATE ( cff_incr_diagnostic )
        ENDIF  ! on STASHflag
        
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',3)

! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,12,STASHwork12,                            &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',4)

        DEALLOCATE (STASHwork12)

      ENDIF !   SF(0,12)

! DEPENDS ON: timer
      If (Ltimer) Call timer ('Diag_adv',4)

      End If ! CycleNo == NumCycles

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS2 S-L Advection',6)

! ----------------------------------------------------------------------
! Section 3.0  Call Atmospheric Physics2
! ----------------------------------------------------------------------
      If (.NOT. L_run_with_physics2) Then
        L_physics_store=l_physics
        L_physics=.false.
      end if

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS3 Atmos_Phys2',5)

      If (ErrorStatus  ==  0 ) Then

        If (L_idealised_data) Then

          If ( CycleNo == 1 ) Then

            ALLOCATE (q_inc_subs(row_length, rows, model_levels))
            ALLOCATE (th_inc_subs(row_length, rows, model_levels))
            ALLOCATE (q_inc_ls(row_length, rows, model_levels))
            ALLOCATE (th_inc_ls(row_length, rows, model_levels))
            ALLOCATE (u_inc_dmp(row_length, rows, model_levels))
            ALLOCATE (q_inc_dmp(row_length, rows, model_levels))
            ALLOCATE (th_inc_dmp(row_length, rows, model_levels))
            ALLOCATE (v_inc_dmp(row_length, n_rows, model_levels))

          End If
! include idealised forcing for theta.

! DEPENDS ON: idl_force
          CALL IDL_Force(                                               &
     &                   R, g, kappa, pi, p_zero                        &
     &,                  row_length, rows, model_levels, timestep       &
     &,                  delta_lambda, delta_phi, model_domain          &
     &,                  lambda_half_width, phi_half_width              &
     &,                  p_max, p_top, p_bottom, p_arbitrary            &
     &,                  lambda_heat_centre, phi_heat_centre            &
     &,                  max_heat_per_day, newtonian_timescale          &
! Dynamical core settings
     &,                  SuHe_newtonian_timescale_ka                    &
     &,                  SuHe_newtonian_timescale_ks                    &
     &,                  SuHe_pole_equ_deltaT, SuHe_static_stab         &
     &,                  SuHe_level_weight, SuHe_sigma_cutoff           &
     &,                  L_SH_Williamson, SuHE_relax                    &
     &,                  L_damp, L_geo_for                              &
     &,                  L_bomex                                        &
     &,                  DMPTIM, HDMP, ZDMP                             &
     &,                  u_geo, v_geo                                   &
     &,                  offx, offy, datastart, at_extremity            &
     &,                  EXNER_THETA_LEVELS                             &
     &,                  P, theta_star                                  &
     &,                  P_THETA_LEVELS, PSTAR                          &
     &,                  THETA                                          &
     &,                  cos_theta_latitude,sin_theta_latitude          &
     &,                  cool_rate, theta_surface                       &
     &,                  timestep_number                                &
     &,                  max_model_levels, max_num_force_times          &
     &,                  Q, q_star                                      &
     &,                  U, V, R_u, R_v                                 &
     &,                  u_ref, v_ref, theta_ref, n_rows                &
     &,                  q_ref                                          &
     &,                  eta_theta_levels, eta_rho_levels               &
     &,                  height_domain                                  &
     &,                  global_row_length, global_rows                 &
     &,                  gc_all_proc_group, nproc                       &
     &,                  tforce_option, qforce_option, uvforce_option   &
     &,                  num_tforce_times, num_qforce_times             &
     &,                  num_uvforce_times                              &
     &,                  tforce_time_interval, qforce_time_interval     &
     &,                  uvforce_time_interval                          &
     &,                  tforce_data_modlev, qforce_data_modlev         &
     &,                  uforce_data_modlev, vforce_data_modlev         &
     &,                  r_rho_levels, r_theta_levels                   &
     &,                  halo_i, halo_j                                 &
     &,                  q_inc_subs, th_inc_subs                        &
     &,                  q_inc_ls, th_inc_ls                            &
     &,                  u_inc_dmp, q_inc_dmp, th_inc_dmp               &
     &,                  v_inc_dmp                                      &
     &,                  f3_at_u, f3_at_v                               &
     &,                  L_physics, problem_number, L_force, LTimer)


          If ( CycleNo == NumCycles ) Then
            DEALLOCATE (q_inc_subs)
            DEALLOCATE (th_inc_subs)
            DEALLOCATE (q_inc_ls)
            DEALLOCATE (th_inc_ls)
            DEALLOCATE (u_inc_dmp)
            DEALLOCATE (q_inc_dmp)
            DEALLOCATE (th_inc_dmp)
            DEALLOCATE (v_inc_dmp)
          End If

        End If      !  L_idealised_data

        If (L_tracer) then

! protect from multiple mem allocations
! save input fields to obtain Atmos_Physics2 increments

! tracers only at final cycle
          If ( CycleNo == NumCycles ) Then
! DEPENDS ON: tr_set_phys
            call TR_Set_Phys(                                           &
                         super_array_size, super_tracer_phys2,          &
                         L_CO2_interactive, CO2,                        &
                         L_Murk_advect, murk,                           &
                         L_Soot, soot_new,                              &
                                 soot_agd,                              &
                                 soot_cld,                              &
                         L_SULPC_SO2, SO2,                              &
                                      SO4_aitken,                       &
                                      so4_accu,                         &
                                      so4_diss,                         &
                         L_sulpc_nh3, nh3,                              &
                         L_sulpc_dms, dms,                              &
                         L_dust, DUST_DIV1,                             &
                                 DUST_DIV2,                             &
                                 DUST_DIV3,                             &
                                 DUST_DIV4,                             &
                                 DUST_DIV5,                             &
                                 DUST_DIV6,                             &
                         L_biomass, bmass_new,                          &
                                    bmass_agd,                          &
                                    bmass_cld,                          &
                         L_ocff, ocff_new,                              &
                                 ocff_agd,                              &
                                 ocff_cld,                              &
                         L_USE_CARIOLLE, OZONE_TRACER,                  &
                         tracer_phys1, tracer, tracer_ukca,             &
                         row_length, rows,                              &
                         model_levels, tr_levels, tr_vars, tr_ukca,     &
                         offx, offy, model_domain,                      &
                         .false., 0, 0                                 &
                         )
          End If ! CycleNo == NumCycles

        end if  ! L_tracer

! save star fields to obtain increments after call to Atmos_Physics2
        If (L_moist_nonhydro_conserve)then
! after SL_thermo q_star holds q_dash (latest estimate to q_(n+1))
          IF(l_mix_ratio)then
! protect from multiple mem allocations
            If ( CycleNo == 1 ) Then
            allocate ( mix_v_phys2(1-offx:row_length+offx,              &
     &                             1-offy:rows+offy, wet_levels) )
            allocate ( mix_cl_phys2(1-offx:row_length+offx,             &
     &                              1-offy:rows+offy, wet_levels) )
            allocate ( mix_cf_phys2(1-offx:row_length+offx,             &
     &                              1-offy:rows+offy, wet_levels) )
            if(L_mcr_qcf2)then
              allocate ( mix_cf2_phys2(1-offx:row_length+offx,          &
     &                                 1-offy:rows+offy, wet_levels) )
            else
              allocate ( mix_cf2_phys2(1,1,1) )
            endif
            if(L_mcr_qrain)then
              allocate ( mix_rain_phys2(1-offx:row_length+offx,         &
     &                                  1-offy:rows+offy, wet_levels) )
            else
              allocate ( mix_rain_phys2(1,1,1) )
            endif
            if(L_mcr_qgraup)then
              allocate ( mix_graup_phys2(1-offx:row_length+offx,        &
     &                                   1-offy:rows+offy, wet_levels) )
            else
              allocate ( mix_graup_phys2(1,1,1) )
            endif
            End If ! CycleNo == 1
! DEPENDS ON: q_to_mix
            call q_to_mix (row_length, rows, wet_levels,                &
     &               offx,offy     ,                                    &
     &               q_star, qcl_star, qcf_star,                        &
     &               qcf2_star, qrain_star, qgraup_star,                &
!    &               .false. ,.false. ,.false.,
     &               L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,             &
     &               mix_v_phys2, mix_cl_phys2, mix_cf_phys2,           &
     &               mix_cf2_phys2, mix_rain_phys2, mix_graup_phys2     &
     &               )
            if(L_pc2)then
! protect from multiple mem allocations
            If ( CycleNo == 1 ) Then
            allocate ( cf_phys2 (1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            allocate ( cfl_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            allocate ( cff_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            End If ! CycleNo == 1
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  cf_phys2(i,j,k)  = cf_star(i,j,k)
                  cfl_phys2(i,j,k) = cfl_star(i,j,k)
                  cff_phys2(i,j,k) = cff_star(i,j,k)
                End Do
              End Do
            End Do
            else
              If ( CycleNo == 1 ) Then
                allocate ( cf_phys2(1,1,1) )
                allocate ( cfl_phys2(1,1,1) )
                allocate ( cff_phys2(1,1,1) )
              End If
            endif
          else

! protect from multiple mem allocations
            If ( CycleNo == 1 ) Then
            allocate ( q_phys2(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
            allocate ( qcl_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            allocate ( qcf_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            if(L_pc2)then
            allocate ( cf_phys2(1-offx:row_length+offx,                 &
     &                          1-offy:rows+offy, wet_levels) )
            allocate ( cfl_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            allocate ( cff_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            else
              allocate ( cf_phys2(1,1,1) )
              allocate ( cfl_phys2(1,1,1) )
              allocate ( cff_phys2(1,1,1) )
            endif
            End If ! CycleNo == 1
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  q_phys2(i,j,k) = q_star(i,j,k)
                  qcl_phys2(i,j,k) = qcl_star(i,j,k)
                  qcf_phys2(i,j,k) = qcf_star(i,j,k)
                End Do
              End Do
            End Do
            if(L_pc2)then
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    cf_phys2(i,j,k)  = cf_star(i,j,k)
                    cfl_phys2(i,j,k) = cfl_star(i,j,k)
                    cff_phys2(i,j,k) = cff_star(i,j,k)
                  End Do
                End Do
              End Do
            endif
           if(L_mcr_qcf2)then
             If ( CycleNo == 1 )                                        &
     &       allocate ( qcf2_phys2(1-offx:row_length+offx,              &
     &                             1-offy:rows+offy, wet_levels) )
             Do k = 1, wet_levels
               Do j = 1, rows
                 Do i = 1, row_length
                   qcf2_phys2(i,j,k)  = qcf2_star(i,j,k)
                 End Do
               End Do
             End Do
           else
             If ( CycleNo == 1 ) allocate ( qcf2_phys2(1,1,1) )
           endif
           if(L_mcr_qrain)then
             If ( CycleNo == 1 )                                        &
     &       allocate ( qrain_phys2(1-offx:row_length+offx,             &
     &                              1-offy:rows+offy, wet_levels) )
             Do k = 1, wet_levels
               Do j = 1, rows
                 Do i = 1, row_length
                   qrain_phys2(i,j,k)  = qrain_star(i,j,k)
                 End Do
               End Do
             End Do
           else
             If ( CycleNo == 1 ) allocate ( qrain_phys2(1,1,1) )
           endif
           if(L_mcr_qgraup)then
             If ( CycleNo == 1 )                                        &
     &       allocate ( qgraup_phys2(1-offx:row_length+offx,            &
     &                               1-offy:rows+offy, wet_levels) )
             Do k = 1, wet_levels
               Do j = 1, rows
                 Do i = 1, row_length
                   qgraup_phys2(i,j,k)  = qgraup_star(i,j,k)
                 End Do
               End Do
             End Do
           else
             If ( CycleNo == 1 ) allocate ( qgraup_phys2(1,1,1) )
           endif
          endif            !L_mix_ratio

        end if  !L_moist_nonhydro_conserve

        If (L_Physics) Then

          If (.NOT. (L_mix_ratio .and. L_moist_nonhydro_conserve)       &
                     .and. CycleNo == 1) Then
            Allocate ( mix_v_phys2 (1-offx:row_length+offx,             &
                         1-offy:rows+offy, wet_levels) )
            Allocate ( mix_cl_phys2(1-offx:row_length+offx,             &
                         1-offy:rows+offy, wet_levels) )
            Allocate ( mix_cf_phys2(1-offx:row_length+offx,             &
                         1-offy:rows+offy, wet_levels) )
            If(L_mcr_qcf2)Then
              Allocate ( mix_cf2_phys2(1-offx:row_length+offx,          &
                                   1-offy:rows+offy, wet_levels) )
            Else
              Allocate ( mix_cf2_phys2(1,1,1) )
            End If
            If(L_mcr_qrain)Then
              Allocate ( mix_rain_phys2(1-offx:row_length+offx,         &
                                   1-offy:rows+offy, wet_levels) )
            Else
              Allocate ( mix_rain_phys2(1,1,1) )
            End If
            If(L_mcr_qgraup)Then
              Allocate ( mix_graup_phys2(1-offx:row_length+offx,        &
                                   1-offy:rows+offy, wet_levels) )
            Else
              Allocate ( mix_graup_phys2(1,1,1) )
            End If

! convert q_star,qcl_star,qcf_star to mix_v_phys2,mix_cl_phys2,mix_cf_phys2
! DEPENDS ON: q_to_mix
            call q_to_mix (row_length, rows, wet_levels,                &
                     offx,offy,                                         &
                     q_star, qcl_star, qcf_star,                        &
                     qcf2_star, qrain_star, qgraup_star,                &
                     L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,             &
                     mix_v_phys2, mix_cl_phys2, mix_cf_phys2,           &
                     mix_cf2_phys2, mix_rain_phys2, mix_graup_phys2     &
                     )
          End If ! (.NOT. (L_mix_ratio .and. L_moist_nonhydro_conserve))

! NB: the star variables and R_u and R_v have not been set in the
!     halo region yet.

! Apply diagnostics only at last cycle.
         If ( CycleNo == NumCycles ) Then

! Allocate diagnostic space for STASH
      IF(L_bl)then
          ALLOCATE (STASHwork3(STASH_maxlen(3,A_im)))
      endif
          IF(convection_option == 3) THEN
            ALLOCATE (STASHwork5(STASH_maxlen(5,A_im)))
          ENDIF   ! on convection_option
      If (L_hydrology)then
          ALLOCATE (STASHwork8(STASH_maxlen(8,A_im)))
      endif
          ALLOCATE (STASHwork9(STASH_maxlen(9,A_im)))
          ALLOCATE (STASHwork19(STASH_maxlen(19,A_im)))
          ALLOCATE (STASHwork26(STASH_maxlen(26,A_im)))

         End If  ! CycleNo == NumCycles

! DEPENDS ON: timer
          If (Ltimer) Call timer ('Atmos_Physics2',3)

! NB if you are changing the argument list to atmos_physics2, please
! do an equivalent change to the two places in routine scm_main to keep
! the single column model consistent.

          IF (nice  ==  1)THEN
! Set sea ice catagory scheme D1 pointers as catogories are not present
! if nice = 1.
            p_ti        => ti
            p_ice_fract => ice_fraction
            p_ice_thick => ice_thickness
          ELSE
            p_ti        => ti_cat
            p_ice_fract => ice_fract_cat
            p_ice_thick => ice_thick_cat
          ENDIF


          If (L_RHCPT) Then
!     Dimension diagnostic 3D RHcrit array
            rhc_row_length = row_length
            rhc_rows = rows
          Else
!     RHcrit will be a 1D parametrized array input from user interface
            rhc_row_length = 1
            rhc_rows = 1
          End if    ! L_RHCPT

         If ( CycleNo == 1 .and. L_physics ) Then
           ALLOCATE ( RHCPT(rhc_row_length, rhc_rows, wet_levels) )
         End If
       If ( CycleNo == 1 ) Then
         If (L_subfilter_horiz .or. L_subfilter_vert) then
           ALLOCATE (FM_3D(row_length,rows,BL_LEVELS))
           ALLOCATE (FH_3D(row_length,rows,BL_LEVELS))
         Else
           ALLOCATE (FM_3D(1,1,1))
           ALLOCATE (FH_3D(1,1,1))
         Endif
       End If

!kdcorbin, 05/10 - initialize mass on first timestep
if (l_co2_mass .or. l_tracer_mass) then
   if (l_tracer_mass .and. i_tracermass_start .lt. 1) then
       i_tracermass_start=1
   endif
   if (l_tracer_mass .and. i_tracermass_start .gt. tr_vars) then
       l_tracer_mass=.false.
       write(6,*) '***********************************'
       write(6,*) 'Tracer Mass Error: '
       write(6,*) '  Incorrect assignment of starting tracer'
       write(6,*) 'Turning off tracer mass fixer'
       write(6,*) 'Specify I_TRACERMASS_START in CNTLATM'
   endif

 if (first_atmstep_call) then
   !DEPENDS ON: TRACER_MASSINIT
   Call TRACER_MASSINIT(     &
  &     row_length, rows, model_levels,   &
  &     tr_levels, tr_vars,               &
  &     halo_i, halo_j, offx, offy, mype, timestep_number,  &
  &     r_theta_levels, r_rho_levels, exner_theta_levels,   &
  &     FV_cos_theta_latitude, delta_lambda, delta_phi,     &
  &     RHO, q, qcl, qcf,                                   &
!CO2 and tracer mass flags
  &     L_CO2_MASS,L_TRACER_MASS,                           &
!CO2 and tracer 3-D fields
  &     CO2, tracer,                                        &
  &     tmass)
   endif
endif

! NB if you are changing the argument list to atmos_physics2, please
! do an equivalent change in routine scm_main to keep the single column
! model consistent.  Note also that there are two calls to
! atmos_physics2 in scm_main.

! DEPENDS ON: atmos_physics2
          Call Atmos_Physics2(halo_i,halo_j,offx,offy,global_row_length &
! Parallel variables
     &, global_rows, gc_proc_row_group, gc_proc_col_group, at_extremity &
     &, nproc, nproc_x, nproc_y, neighbour, g_rows, g_row_length        &
     &,     g_datastart, mype, NumCycles, CycleNo                       &

! field dimensions etc.
     &,     row_length, rows, n_rows, land_field, model_levels, nice    &
     &,     wet_levels, bl_levels, st_levels, sm_levels, cloud_levels   &
     &,     land_ice_points, soil_points, n_cca_lev, ntiles, tr_levels  &
     &,     first_constant_r_rho_level, DIM_CS1, DIM_CS2                &

! IN Substepping information and switches
     &,     Num_Substeps, L_phys2_substep, L_regular                    &
     &,     l_mr_physics2, model_domain, L_dry                          &
     &,     FORMDRAG, OROG_DRAG_PARAM, LCAL360, L_emcorr, Ltimer        &
     &,     L_us_blsol,BL_OPTIONS,L_CLD_AREA,L_ACF_Cusack,L_ACF_Brooks  &
     &,     L_RHCPT,L_hydrology,L_bl,L_Murk,L_murk_advect,L_murk_source &
     &,     L_BL_TRACER_MIX,L_DUST,L_CAM_DUST,L_SULPC_SO2,L_SULPC_NH3   &
     &,     L_SBLeq,L_SBLco,L_sulpc_dms,L_soot,L_biomass,L_ocff         &
     &,     L_co2_interactive,L_ctile, L_co2_emits, L_use_bl_diag_term  &
     &,     L_pc2,L_pc2_reset,L_eacf,L_LAMBDAM2,L_FULL_LAMBDAS,L_ukca   &
     &,     L_sice_heatflux, L_OHADGEM1, L_USE_CARIOLLE                 &
     &,     l_anthrop_heat_src                                          &
! Model Parameters
     &,     alpha_Cd, Puns, Pstb, RHcrit, CO2_MMR, tr_vars, tr_ukca     &
     &,     Muw_SBL, Mwt_SBL, cloud_fraction_method, overlap_ice_liquid &
     &,     ice_fraction_method, ctt_weight, t_weight, qsat_fixed       &
     &,     sub_cld, x1i, x1ic, x1r, x2r, x4r, l_psd, ai, bi, aic, bic  &
     &,     lsp_ei,lsp_fi,lsp_eic,lsp_fic                               &
     &,     dbsdtbs_turb_0, Charnock, SeaSalinityFactor                 &
! Physical constants
     &,Lc,Lf,Cp, two_Omega,p_zero, kappa, R, g, Lapse, earth_radius, Pi &
! Vertical coordinate levels.
     &,r_rho_levels, r_theta_levels, r_at_u, r_at_v, unscaled_dry_rho   &
     &,     eta_theta_levels, eta_rho_levels, delta_lambda, delta_phi   &
     &,     gdlambda_p(lambda_start), dphi_p                            &
     &,     wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v                &
     &,     lat_rot_NP, long_rot_NP, f3_at_u                            &
! Variables required by STPH_SCV
! add mix ratios 
     &,     mix_v_phys2, mix_cl_phys2, mix_cf_phys2                     &
! Time stepping information and trig arrays
     &, timestep, I_year, I_day_number, I_hour, I_minute, I_second      &
     &, timestep_number, sin_theta_longitude, cos_theta_longitude       &
!     &,     FV_cos_theta_latitude,sec_theta_latitude,                   &
     &, sin_theta_latitude,FV_cos_theta_latitude,sec_theta_latitude,    &
! River routing
     &AOCPL_ROW_LENGTH,AOCPL_P_ROWS,L_RIVERS,XPA,XUA,XVA,YPA,YUA,YVA,   &
     &G_P_FIELD,G_R_FIELD,A_INTHD(16),lasize(1,fld_type_r,              &
     &halo_type_no_halo),lasize(2,fld_type_r,halo_type_no_halo),        &
     &glsize(1,fld_type_r),glsize(2,fld_type_r),RIVER_STEP,RIVER_VEL,   &
     &RIVER_MCOEF,RIV_DIRECTION,RIV_SEQUENCE,RIV_STORAGE,               &
!  Add inland basin outflow to arguments
     &  RIV_INLANDATM,L_INLAND,                                         &


! Grid-to-grid river routing
     &  RIV_IAREA,RIV_SLOPE, RIV_FLOWOBS1,RIV_INEXT                     &
     &  ,RIV_JNEXT,RIV_LAND,RIV_SUBSTORE,                               &
     &  RIV_SURFSTORE,RIV_FLOWIN,RIV_BFLOWIN,                           &

!
! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &      STASHwork3,STASHwork5,STASHwork8,STASHwork9,STASHwork19,    &
     &      STASHwork26,                                                &
!
! SCM Diagnostics (dummy values in full UM)
     &  nSCMDpkgs, L_SCMDiags,                                          &
!
! Data Fields.
     & THETA, Q, QCL, QCF,RHO                                           &
     &,U, V, W, W_ADV,P                                                 &

     &,PSTAR, EXNER_RHO_LEVELS,EXNER_THETA_LEVELS                       &
     &,LAND, P_THETA_LEVELS                                             &
! variables for subgrid turbulence scheme
     &,visc_BL_m, FM_3D, FH_3D, L_subfilter_vert, L_subfilter_horiz     &
     &,L_subfilter_blend,max_diff, turb_startlev_vert,turb_endlev_vert  &
     &,BL_COEF_KM,BL_COEF_KH                                            &

! ancillary fields and fields needed to be kept from timestep to
! timestep

     &,land_index, land_ice_index, soil_index, CANOPY_WATER             &
     &,SNODEP, THERM_COND, THERM_CAP, VOL_SMC_CRIT                      &
     &,VOL_SMC_WILT, VOL_SMC_SAT, STHF, STHU                            &
     &,OROG_SIL,OROG_HO2,ICE_THICKNESS,ICE_FRACTION                     &
     &,U_SEA,V_SEA,U_0_P,V_0_P,CCA,CCB                                  &
     &,CCT, CCLWP, CCW_RAD, LCBASE, DEEP_SOIL_TEMP, p_ti, TI            &
     &,TSTAR,Z0, p_ice_fract, p_ice_thick                               &
     &,SAT_SOIL_COND,SAT_SOILW_SUCTION,CLAPP_HORN                       &
     &,SMCL, T1_SD, Q1_SD, ZH                                           &
     &,CF_AREA, CF_BULK, CF_LIQUID                                      &
     &,CF_FROZEN, ls_rain, ls_snow, micro_tends                         &
     &,photosynth_act_rad, rad_hr, surf_radflux, SOIL_CLAY              &
     &,SOIL_SILT,SOIL_SAND, DUST_MREL1,DUST_MREL2                       &
     &,DUST_MREL3,DUST_MREL4,DUST_MREL5,DUST_MREL6                      &
     &,SO2_HILEM, SO2_EM, NH3_EM, DMS_EM                                &
     &,SOOT_HILEM, SOOT_EM, OCFF_HILEM, OCFF_EM, CO2_EMITS, CO2FLUX     &

! tracer fluxes - kdcorbin, 05/10
     &, TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4          &
     &, TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8          &
     &, TRACER_FLUX9, TRACER_FLUX10,TRACER_FLUX11,TRACER_FLUX12         &
     &, TRACER_FLUX13,TRACER_FLUX14,TRACER_FLUX15,TRACER_FLUX16         &
     &, TRACER_FLUX17,TRACER_FLUX18,TRACER_FLUX19,TRACER_FLUX20         &
! CO2 global emissions for mass balance - kdcorbin, 05/10               &
     &, CO2EMITMASS                                                     &
! rml 1/7/13 flag for co2 flux into passive tracer
     &, L_CO2_TRACER                                                    &

! IN/OUT
     &,theta_star,q_star,qcl_star,qcf_star,cf_star,cfl_star,cff_star    &
     &,R_u, R_v, R_w                                                    &
     &,net_flux,net_mflux,murk,tracer,tracer_ukca                       &
     &, DUST_DIV1, DUST_DIV2, DUST_DIV3                                 &
     &, DUST_DIV4, DUST_DIV5, DUST_DIV6                                 &
     &, so2, dms, so4_aitken, so4_accu                                  &
     &, so4_diss, nh3, soot_new                                         &
     &, soot_agd, soot_cld, bmass_new                                   &
     &, bmass_agd, bmass_cld, ocff_new, ocff_agd, ocff_cld, co2         &
! IN/OUT STPH_RP
     &,G0_RP,par_mezcla                                                 &

! IN/OUT River routing
     &, TOT_SURFROFF, TOT_SUBROFF                                       &
!
! OUT Fields
     &, rho_km, cH, NTML, CUMULUS, NBDSC, NTDSC                         &
     &, rhcpt, rhc_row_length, rhc_rows                                 &

! Additional variables for MOSES II
     &, FRAC_TYP, DISTURB_VEG, CANHT_PFT, LAI_PFT                       &
     &, CAN_WATER_TILE, CATCH_TILE, CATCH_SNOW                          &
     &, SNOW_GRND, SNODEP_TILE, Z0_TILE, TSTAR_TILE                     &
     &, INFIL_TILE, RGRAIN_TILE, CS, GS                                 &
     &, co2_dim_row, co2_dim_len, l_neg_tstar, l_snow_albedo            &
     &, l_phenol, l_triffid, l_trif_eq, l_q10, A_INTHD(23)              &
     &, STEPim(atmos_im),phenol_period,A_INTHD(22),CAN_MODEL            &
     &, G_LF_PFT_ACC, G_PHLF_PFT_ACC, NPP_PFT_ACC                       &
     &, RSP_W_PFT_ACC, RSA                                              &
     &, land_pts_trif, npft_trif, dolr, LW_down, SW_tile                &
     &,FRAC_LAND,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE                        &
     &,SOIL_ALB,cos_zenith_angle,can_rad_mod,ilayers                    &
     &, RADNET_TILE                                                     &
!glr cable start-----------------
     &, surf_down_sw,alb_tile,l_tile_pts                                &
     &, lat,cos_theta_longitude,day,time_sec,SW_DOWN                    &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,SNOW_RHO3L        &
     &, SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,TSOIL_TILE,T_SURF_TILE &
     &, HCONS,SOIL_TYPE,VEG_TYPE                                        &
     &, SNOW_FLG3L,total_nsteps                                         &
     &           ,FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB               &
     &           ,TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB           &
     &           ,USTAR_CAB,SURF_HTF_CAB                                &
! Lestevens - need first_atmstep_call for wblake fix
     &, l_cable,wblake_ratio,WB_LAKE                                    &
     &, TOT_WBLAKE,TOT_SUBRUN                                           &
     &           ,TOT_ALB                                               &
     &           ,U_S_CAB,CH_CAB,CD_CAB                                 &
     &           ,CD,CHX                                                &
     &           ,TILE_PTS,TILE_INDEX                                   &  
     &, SNOW_AGE,RTSOIL_TILE                                            &
     &, GFLUX_TILE,SGFLUX_TILE                                          &
     &, TRANSP_TILE                                                     &
     &, CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER                     &
     &, NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE                           &
!glr cable end -----------------
 
! Additional variables required for large-scale hydrology:
     &, L_TOP,L_PDM,FEXP,GAMMA_INT,TI_MEAN,TI_SIG                       &
     &, FSFC_SAT,F_WETLAND,WATER_TABLE,STHZW                            &
     &, A_FSAT,C_FSAT,A_FWET,C_FWET,L_SOIL_SAT_DOWN                     &
! Cariolle ozone 
     &, OZONE_TRACER                                                    &

! SCM and idealised UM surface forcing parameters
     &,L_flux_bc,flux_e,flux_h,L_spec_z0,z0m_scm,z0h_scm                &
!
! error information
     &,                      ErrorStatus,                               &
      !end step of experiment
     endstep, mype)

! DEPENDS ON: timer
          If (Ltimer) Call timer ('Atmos_Physics2',4)

          If(.NOT. (L_mix_ratio .and. L_moist_nonhydro_conserve)        &
     &            .and. CycleNo == NumCycles)then
            Deallocate(mix_graup_phys2)
            Deallocate(mix_rain_phys2)
            Deallocate(mix_cf2_phys2)
            Deallocate(mix_cf_phys2)
            Deallocate(mix_cl_phys2)
            Deallocate(mix_v_phys2)
          End If

!-----------------------------
!Tracer/CO2 Mass Fixer
!kdcorbin, 05/10
!-----------------------------

if (l_co2_mass .or. l_tracer_mass) then
! DEPENDS ON: tracer_massfix
  call tracer_massfix(row_length,rows,model_levels                    &
             ,tr_levels,tr_vars,i_tracermass_start                    &
             ,halo_i, halo_j, offx, offy, mype                        &
             ,timestep,timestep_number                                &
             ,r_theta_levels, r_rho_levels, exner_theta_levels        &
             ,FV_cos_theta_latitude, delta_lambda, delta_phi          &
             ,RHO, q, qcl, qcf                                        &
             ,L_CO2_MASS,L_TRACER_MASS,CO2, tracer               &
             ,co2emitmass,tracer_flux1,tracer_flux2,tracer_flux3        &
             ,tracer_flux4,tracer_flux5,tracer_flux6,tracer_flux7     &
             ,tracer_flux8,tracer_flux9,tracer_flux10,tracer_flux11   &
             ,tracer_flux12,tracer_flux13,tracer_flux14,tracer_flux15 &
             ,tracer_flux16,tracer_flux17,tracer_flux18,tracer_flux19 &
             ,tracer_flux20,tmass)

endif 

!-----------------------
!Methane Loss via OH
!kdcorbin, 05/10
!-----------------------

if (tr_vars > 0) then
!if (l_methane_loss .or. l_mcf_loss) then
   allocate(lossrate(row_length,rows,tr_levels,tr_vars))
   lossrate=0.
endif

if (l_methane_loss) then

  if (mype .eq. 0) then
     write(6,*) 'Calculating Methane Loss via OH, O1D and CL'
  endif

  if ((i_methane_tracers .lt. 1) .or.  &
     (i_methane_tracers .gt. tr_vars)) then
      write(6,*) '**********************************'
      write(6,*) 'Methane Specification Error: '
      write(6,*) '  Incorrect assignment of methane tracers'
      write(6,*) 'Specify I_METHANE_TRACERS in CNTLATM'   
      write(6,*) 'Stopping model'
      stop
   endif

! DEPENDS ON: tracer_methaneloss
  call tracer_methaneloss(row_length,rows,model_levels,i_methane_tracers  &
             ,tr_levels,tr_vars,offx,offy,timestep       &
             ,theta,p,rho,tracer,oh,ho2,lossrate)
endif  !methane loss

!-----------------------
!Radon Decay
!kdcorbin, 05/10
!-----------------------
if (l_radon_decay) then
 
   if (mype .eq. 0) then
      write(6,*) 'Calculating Radon Decay'
   endif

  if ((i_radon_tracernumber .lt. 1) .or.  &
      (i_radon_tracernumber .gt. tr_vars)) then
      write(6,*) '**********************************'
      write(6,*) 'Radon Specification Error: '
      write(6,*) '  Incorrect assignment of radon tracer'
      write(6,*) 'Specify I_RADON_TRACERNUMBER in CNTLATM'   
      write(6,*) 'Stopping model'
      stop
   endif

! DEPENDS ON: tracer_radondecay
  call tracer_radondecay(row_length,rows,tr_levels,tr_vars  &
          ,i_radon_tracernumber,offx,offy,timestep,tracer)

endif

!-----------------------
!MCF Loss
!kdcorbin, 05/10
!-----------------------
if (l_mcf_loss) then

   if (mype .eq. 0) then
      write(6,*) 'Calculating MCF Loss'
   endif

  if ((i_mcf_tracernumber .lt. 1) .or.  &
      (i_mcf_tracernumber .gt. tr_vars)) then
      write(6,*) '**********************************'
      write(6,*) 'MCF Specification Error: '
      write(6,*) '  Incorrect assignment of MCF tracer'
      write(6,*) 'Specify I_MCF_TRACERNUMBER in CNTLATM'   
      write(6,*) 'Stopping model'
      stop
   endif

! DEPENDS ON: tracer_mcfloss
  call tracer_mcfloss(row_length,rows,model_levels   &
         ,i_mcf_tracernumber   &
         ,tr_levels,tr_vars,offx,offy,timestep       &
         ,halo_i,halo_j,r_theta_levels,r_rho_levels  &
         ,theta,p,rho,tracer,oh,h2o2_limit,o3_chem,lossrate)

endif

!------------------------
!Tracer Information Printout
!kdcorbin, 05/10
!-----------------------
if (tr_vars > 0) then
!if (l_methane_loss .or. l_mcf_loss) then
   !DEPENDS ON: tracer_massprint
   call tracer_massprint(row_length,rows,model_levels                    &
             ,tr_levels,tr_vars                                       &
             ,halo_i, halo_j, offx, offy, mype                        &
             ,timestep, timestep_number                               &
             ,r_theta_levels, r_rho_levels, exner_theta_levels        &
             ,FV_cos_theta_latitude, delta_lambda, delta_phi          &
             ,RHO, q, qcl, qcf, tracer                                &
             ,lossrate,tracer_flux1,tracer_flux2,tracer_flux3        &
             ,tracer_flux4,tracer_flux5,tracer_flux6,tracer_flux7     &
             ,tracer_flux8,tracer_flux9,tracer_flux10,tracer_flux11   &
             ,tracer_flux12,tracer_flux13,tracer_flux14,tracer_flux15 &
             ,tracer_flux16,tracer_flux17,tracer_flux18,tracer_flux19 &
             ,tracer_flux20)

   deallocate(lossrate)
endif


! ----------------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after  Atmos_Physics2 ***'
        ENDIF ! L_print_pe .or. mype ==0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      endif !  L_print_L2norms
! ---------------------------------------------------------------

! Apply diagnostics only at last cycle.
         If ( CycleNo == NumCycles ) Then

! Diagnostics STASHed for each section in Atmos_Physics2:

! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
      IF(L_bl)then
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,3,STASHwork3,                            &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork3)
      endif

          IF(convection_option == 3) THEN
! DEPENDS ON: timer
            IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
            CALL STASH(a_sm,a_im,5,STASHwork5,                          &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
            IF (Ltimer) CALL TIMER('STASH',4)
            DEALLOCATE (STASHwork5)
          ENDIF   ! on convection_option

      If (L_hydrology)then
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,8,STASHwork8,                            &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork8)
      endif
!
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,9,STASHwork9,                            &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork9)
!
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,19,STASHwork19,                          &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork19)
!
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,26,STASHwork26,                          &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork26)
!
          End If  ! CycleNo == NumCycles

          If (.NOT. L_run_with_physics2 ) Then
            L_physics=L_physics_store
          End If

!L --------- UM Section 18---- Data Assimilation --------------------
! DEPENDS ON: timer
          IF (LTimer) CALL TIMER('AS18 Assimilation',5)

          IF(L_AC .AND. LASSIMILATION .AND. ErrorStatus == 0)THEN
! Do AC assimilation

! copy non halo values of _star moisture variables into WORK arrays
! for passing to AC_CTL where they can be interpreted more
! conveniently as 2-d arrays (and halos not required)
            ALLOCATE ( STASHwork18(STASH_maxlen(18,A_im)) )
            ALLOCATE ( work_q(row_length,rows,wet_levels) )
            ALLOCATE ( work_qcl(row_length,rows,wet_levels) )
            ALLOCATE ( work_qcf(row_length,rows,wet_levels) )
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  work_q(i,j,k)   = q_star(i,j,k)
                  work_qcl(i,j,k) = qcl_star(i,j,k)
                  work_qcf(i,j,k) = qcf_star(i,j,k)
                End Do
              End Do
            End Do

! DEPENDS ON: timer
            IF(LTIMER)  CALL TIMER('AC_CTL',3)

! DEPENDS ON: ac_ctl
            CALL AC_CTL(STASH_MAXLEN(18,atmos_im),theta_field_size,     &
     &                  wet_levels, model_levels,                       &
     &                  theta_star, work_q, work_qcl, work_qcf,         &
     & OBS_FLAG,OBS,obs_flag_len,obs_len,                               &
     & p, p_theta_levels, exner_theta_levels,                           &
     & r_theta_levels, FV_cos_theta_latitude,                           &
     & cf_area, cf_bulk, cf_liquid, cf_frozen,                          &
     & pstar, ntml, cumulus,                                            &
     & STASHwork18,                                                     &
! ARGDUMA Dump headers
     &  A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
     &  A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
     &  A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! #include "arglndm.h" 
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &                  l_mr_acctl,ErrorStatus,CMessage)
! DEPENDS ON: timer
            IF(LTIMER) CALL TIMER('AC_CTL',4)

! DEPENDS ON: timer
            IF(LTIMER) CALL TIMER('STASH',3)

! moved from out of AC_CTL
! DEPENDS ON: stash
        CALL STASH(a_sm, a_im, 18, STASHWORK18,                         &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &                               ErrorStatus,CMessage)

! DEPENDS ON: timer
            IF(LTIMER) CALL TIMER('STASH',4)

! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

! copy back new values of _star moisture variables
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  q_star  (i,j,k) = work_q(i,j,k)
                  qcl_star(i,j,k) = work_qcl(i,j,k)
                  qcf_star(i,j,k) = work_qcf(i,j,k)
                End Do
              End Do
            End Do
            DEALLOCATE (work_q)
            DEALLOCATE (work_qcl)
            DEALLOCATE (work_qcf)
            DEALLOCATE (STASHwork18)
          END IF      ! LASSIMILATION and L_AC

! DEPENDS ON: timer
          IF (LTimer) CALL TIMER('AS18 Assimilation',6)

          If (.not. L_use_bl_diag_term) Then
! zero ch term
            Do k = 1, model_levels-1
              Do j = 1-offy, rows+offy
                Do i = 1-offx, row_length+offx
                  cH(i,j,k) = 0.
                End Do
              End Do
            End Do
            Do k = 0, bl_levels-1
              Do j = 1-offy, rows+offy
                Do i = 1-offx, row_length+offx
                  rho_km(i,j,k) = 0.
                End Do
              End Do
            End Do
          End If    ! L_use_bl_diag_term

       Else    ! L_physics .false.

! zero ch term
          Do k = 1, model_levels-1
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                cH(i,j,k) = 0.
              End Do
            End Do
          End Do
! zero rho_km term
          Do k = 0, bl_levels-1
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                rho_km(i,j,k) = 0.
              End Do
            End Do
          End Do

! Include friction terms

! Run with a positive timestep if integrating backwards.
          IF (L_Backwards) timestep = pos_timestep
!  When no physics applied then may need a simple friction
!   if non-inviscid (idealised_problem) then DO NOT apply friction

          if(problem_number  /=  idealised_problem) Then
!  standard simple friction being used
            if(problem_number  ==  dynamical_core) Then
! DEPENDS ON: idl_friction_suarez_held
              Call IDL_Friction_Suarez_Held(                            &
     &                         row_length, rows, n_rows                 &
     &,                        model_levels, timestep                   &
     &,                        model_domain                             &
     &,                        offx, offy, at_extremity                 &
     &,                        friction_level                           &
     &,                        base_frictional_timescale                &
     &,                        SuHe_sigma_cutoff, SuHe_fric             &
     &,                        P, PSTAR                                 &
     &,                        U, V                                     &
     &,                        R_u, R_v )
            Else     ! problem_number  /=  dynamical_core
              Do k = 1, bl_levels
                Do j = j_begin, j_end
                  Do i = 1, row_length
                    R_u(i,j,k) = R_u(i,j,k) - timestep *                &
     &                     friction_level(k) * u(i,j,k)
                  End Do
                End Do
                Do j = 1, n_rows
                  Do i = 1, row_length
                    R_v(i,j,k) = R_v(i,j,k) - timestep *                &
     &                     friction_level(k) * v(i,j,k)
                  End Do
                End Do
              End Do

            EndIf    ! problem_number  ==  dynamical_core

          EndIf    ! problem_number  /=  idealised_problem

! Go back to negative timestep if integrating backwards.
          IF (L_Backwards) timestep = neg_timestep

        End If !   ! L_physics

! remove any instability
        If (L_adjust_wet .or.                                           &
     &      L_dry .or.                                                  &
     &      ( .not. L_physics .and.                                     &
     &        (.not. L_idealised_data .or. problem_number  ==  1        &
     &        .or. problem_number  ==  2) ) )Then

! DEPENDS ON: timer
          If (Ltimer) Call timer ('dry_static_adj',3)
! DEPENDS ON: dry_static_adj
          Call dry_static_adj(                                          &
     &                       theta_star, q_star, epsilon,               &
     &                       rows, row_length, model_levels,            &
     &                       wet_levels, offx, offy,                    &
     &                       Instability_diagnostics,                   &
     &                       .false.)
! DEPENDS ON: timer
          If (Ltimer) Call timer ('dry_static_adj',4)

        End If   ! L_adjust_wet .or. L_dry ....

      End If   ! ErrorStatus = 0

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS3 Atmos_Phys2',6)

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS3 Diffusion',5)

      If (ErrorStatus  ==  0 ) Then

        If ( L_diff_ctl ) then

! Apply diagnostics only at last cycle.
         If ( CycleNo == NumCycles ) Then

! Save current values to form diagnostics increments over diffusion
          IF( sf(185,13)) THEN
            ALLOCATE ( u_incr_diagnostic(row_length,rows,model_levels) )
! Hold u increment (no halos needed)
            DO k=1,model_levels
              DO j=1,rows
                DO i=1,row_length
                  u_incr_diagnostic(i,j,k) =R_u(i,j,k)
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDIF                   ! on STASHflag

          IF( sf(186,13) ) THEN
            ALLOCATE ( v_incr_diagnostic(row_length,n_rows,             &
     &                                              model_levels) )
! Hold v increment (no halos needed)
            DO k=1,model_levels
              DO j=1,n_rows
                DO i=1,row_length
                  v_incr_diagnostic(i,j,k) =R_v(i,j,k)
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDIF                    ! on STASHflag

          IF( sf(181,13) ) THEN
            ALLOCATE ( T_incr_diagnostic(row_length,rows,model_levels) )
! note at this point theta_star holds all increments to theta
            DO k=1,model_levels
              DO j=1,rows
                DO i=1,row_length
                  T_incr_diagnostic(i,j,k) =theta_star(i,j,k)
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDIF                    ! on STASHflag

          IF( sf(182,13) ) THEN
            ALLOCATE ( q_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point q_star holds all increments to q
            DO k=1,wet_levels
              DO j=1,rows
                DO i=1,row_length
                  q_incr_diagnostic(i,j,k) =q_star(i,j,k)
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDIF                    ! on STASHflag

         End If  ! CycleNo == NumCycles
          IF ( CycleNo == 1 ) Then
          IF( sf(201,13) ) THEN
            ALLOCATE ( w_local_mask(row_length,rows) )
          ENDIF                    ! on STASHflag

          Endif  ! CycleNo == 1

          If (L_subfilter_horiz) then
!
! visc_m is currently S * lengthscale^2. Now multiply this by
! stability fn. Both FH and FM are on w points at this point.
!
            Do k = 1, model_levels - 2

              If (k >= turb_startlev_horiz .AND.                        &
     &               k <= turb_endlev_horiz) Then

                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = visc_m(i,j,k)*FH_3D(i,j,k+1)
                    visc_m(i,j,k) = visc_m(i,j,k)*FM_3D(i,j,k+1)
                    visc_h(i,j,k) = min(visc_h(i,j,k), max_diff)
                    visc_m(i,j,k) = min(visc_m(i,j,k), max_diff)
                  End Do
                End Do

              Else

                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = 0.0
                    visc_m(i,j,k) = 0.0
                  End Do
                End Do

              End If

            End Do
!
! Set visc_m and visc_h at the top two levels
!
            Do k = model_levels - 1, model_levels

              If (turb_startlev_horiz <= k .AND.                        &
     &            turb_endlev_horiz >= k) Then
     
                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = visc_h(i,j,model_levels-2)
                    visc_m(i,j,k) = visc_m(i,j,model_levels-2)
                  End Do
                End Do

              Else

                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = 0.0
                    visc_m(i,j,k) = 0.0
                  End Do
                End Do

              End If

            End Do

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(visc_m,                                    &
     &                       row_length, rows, model_levels,            &
     &                       halo_i, halo_j, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(visc_h,                                    &
     &                       row_length, rows, model_levels,            &
     &                       halo_i, halo_j, fld_type_p, .false.)

            If (model_domain  ==  mt_lam .or. model_domain  ==          &
     &                                          mt_cyclic_lam ) Then

! DEPENDS ON: fill_external_halos
              Call Fill_external_halos (visc_m,                         &
     &                                  row_length, rows, model_levels, &
     &                                  halo_i, halo_j)

! DEPENDS ON: fill_external_halos
              Call Fill_external_halos (visc_h,                         &
     &                                  row_length, rows, model_levels, &
     &                                  halo_i, halo_j)

            End If

          Else If (L_subfilter_vert) then

! visc_m is currently S * lengthscale^2. Now multiply this by
! stability fn. Both FH and FM are on w points at this point.
!
            Do k = 1, model_levels - 2

              If (k >= turb_startlev_vert .AND.                        &
     &               k <= turb_endlev_vert) Then
! use stability functions from BL scheme

                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = visc_m(i,j,k)*FH_3D(i,j,k+1)
                    visc_m(i,j,k) = visc_m(i,j,k)*FM_3D(i,j,k+1)
                  End Do
                End Do

              Else

                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = 0.0
                    visc_m(i,j,k) = 0.0
                  End Do
                End Do

              End If

            End Do
!
! Set visc_m and visc_h at the top two levels
!
            Do k = model_levels - 1, model_levels
              If (turb_startlev_vert <= k .AND.                         &
     &            turb_endlev_vert >= k) Then 
                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = visc_h(i,j,model_levels-2)
                    visc_m(i,j,k) = visc_m(i,j,model_levels-2)
                  End Do
                End Do
              Else
                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = 0.0
                    visc_m(i,j,k) = 0.0
                  End Do
                End Do
              End If
            End Do


! DEPENDS ON: swap_bounds
            Call Swap_Bounds(visc_m,                                    &
     &                       row_length, rows, model_levels,            &
     &                       halo_i, halo_j, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(visc_h,                                    &
     &                       row_length, rows, model_levels,            &
     &                       halo_i, halo_j, fld_type_p, .false.)

            If (model_domain == mt_lam .or. model_domain ==             &
     &                                        mt_cyclic_lam ) Then     

! DEPENDS ON: fill_external_halos
              Call Fill_external_halos (visc_m,                         &
     &                                  row_length, rows, model_levels, &
     &                                  halo_i, halo_j)

! DEPENDS ON: fill_external_halos
              Call Fill_external_halos (visc_h,                         &
     &                                  row_length, rows, model_levels, &
     &                                  halo_i, halo_j)

            End If


          End If  ! L_subfilter_horiz or L_subfilter_vert

! DEPENDS ON: ni_diff_ctl
          Call NI_Diff_Ctl(                                             &
     &                     L_diffusion, L_cdiffusion, L_subfilter_horiz,&
     &                     L_vertical_diffusion, L_divdamp,             &
     &                     L_ramp, ramp_lat_radians,                    &
     &                     L_Backwards, Ltimer,                         &
     &                     timestep, pos_timestep, neg_timestep,        &
     &                     THETA, W, Q,                                 &
     &                     QCL, QCF,                                    &
     &                     QCF2, QRAIN, QGRAUP,                         &
     &                     mix_v, mix_cl, mix_cf,                       &
     &                     mix_cf2, mix_rain, mix_graup,                &
     &                     L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,       &
     &                     U, V, RHO,                                   &
     &                     EXNER_RHO_LEVELS, W_ADV,                     &
     &                     r_theta_levels, r_rho_levels, r_at_u, r_at_v,&
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     sec_theta_latitude,                          &
     &                     cos_theta_latitude, sec_v_latitude,          &
     &                     FV_sec_theta_latitude, cos_v_latitude,       &
     &                     sin_theta_latitude, sin_v_latitude, pi,      &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity, gc_proc_row_group,             &
     &                     mype, nproc, nproc_x, nproc_y, neighbour,    &
     &                     delta_lambda, delta_phi, L_regular,          &
     &                     glambda_p(lambda_start), phi_p,              &
     &                     glambda_u(lambda_start), phi_v,              &
     &                     rows, row_length, n_rows,                    &
     &                     model_levels, wet_levels, model_domain,      &
     &                     global_row_length, global_rows,              &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_coefficient_w,                     &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_coefficient_wind,                  &
     &                     diffusion_order_thermo,                      &
     &                     diffusion_order_w,                           &
     &                     diffusion_order_q,                           &
     &                     diffusion_order_wind,                        &
     &                     visc_m, visc_h,                              &
     &                     horizontal_level, tar_horizontal,            &
     &                     level_start_wind, level_stop_wind,           &
     &                     level_start_q, level_stop_q,                 &
     &                     level_start_theta, level_stop_theta,         &
     &                     L_tardiff_q, w_conv_limit, tardiffq_factor,  &
     &                     tardiffq_test, tardiffq_start, tardiffq_end, &
     &                     sf(201,13), w_local_mask,                    &
     &                     L_adjust_theta,                              &
     &                     adjust_theta_start, adjust_theta_end,        &
     &                     L_vdiff_uv, vdiffuv_test, epsilon,           &
     &                     vdiffuv_factor, vdiffuv_start, vdiffuv_end,  &
     &                     vert_diffusion_coeff_wind,                   &
     &                     vert_diffusion_coeff_q,                      &
     &                     vert_diffusion_coeff_theta,                  &
     &                     div_damp_coefficient,                        &
     &                     theta_star, R_w, q_star,                     &
     &                     qcl_star, qcf_star,                          &
     &                     qcf2_star, qrain_star, qgraup_star,          &
     &                     R_u, R_v,                                    &
     &                     L_mix_ratio)

! ----------------------------------------------------------------------
! Section 13.1 Diagnostics at from diffusion and divergence damping
! ----------------------------------------------------------------------
! Apply diagnostics only at last cycle.
        If ( CycleNo == NumCycles ) Then

! DEPENDS ON: timer
          If (Ltimer) Call timer ('Diag_dif',3)

! section 13:
          IF( SF(0,13)                                                  &
                               ! Diagnostics required for this section
     &        .AND. ErrorStatus == 0) THEN

! Allocate diagnostic space for STASH
            IF(.NOT. L_Filter)then
            ALLOCATE (STASHwork13(STASH_maxlen(13,A_im)))
            ENDIF

! DEPENDS ON: diagnostics_dif
            CALL Diagnostics_dif(                                       &
     &         row_length, rows, n_rows, model_levels, wet_levels,      &
! primary  fields:
     &         THETA, Q,                                                &
! wind field increments after  dif :
     &         R_u, R_v,                                                &
! wind field increments before diffusion (on stashflag):
     &         u_incr_diagnostic, v_incr_diagnostic,                    &
     & T_incr_diagnostic,q_incr_diagnostic,                             &
! variables for subgrid turbulence scheme
     &         L_subfilter_horiz, L_subfilter_vert, visc_m, visc_h,     &
     &         shear, RNEUTML, FM_3D, FH_3D,                            &
     &         BL_LEVELS,BL_COEF_KM,BL_COEF_KH,                         &
     &         w_local_mask,                                            &
     &         theta_star, q_star,                                      &
     &         EXNER_THETA_LEVELS,                                      &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &         STASHwork13)

! Tidy allocatable arrays
            IF( sf(185,13) ) THEN
              DEALLOCATE ( u_incr_diagnostic )
            ENDIF  ! on STASHflag

            IF( sf(201,13) ) THEN
              DEALLOCATE ( w_local_mask )
            ENDIF  ! on STASHflag

            IF( sf(186,13) ) THEN
              DEALLOCATE ( v_incr_diagnostic )
            ENDIF  ! on STASHflag

            IF( sf(181,13) ) THEN
              DEALLOCATE ( T_incr_diagnostic )
            ENDIF  ! on STASHflag

            IF( sf(182,13) ) THEN
              DEALLOCATE ( q_incr_diagnostic )
            ENDIF  ! on STASHflag

            IF( .NOT. ( (Model_domain==mt_global ) .AND.               &
     &                (L_polar_filter_incs .or. L_filter_incs) ) )THEN

! DEPENDS ON: timer
            IF (Ltimer) CALL TIMER('STASH',3)

! DEPENDS ON: stash
            CALL STASH(a_sm,a_im,13,STASHwork13,                        &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &      ErrorStatus,Cmessage)

! DEPENDS ON: timer
            IF (Ltimer) CALL TIMER('STASH',4)

            DEALLOCATE (STASHwork13)
            
          ENDIF ! .NOT. ( (Model_domain==mt_global ) etc

          ENDIF !   SF(0,13)

! DEPENDS ON: timer
          If (Ltimer) Call timer ('Diag_dif',4)

          End If ! CycleNo == NumCycles
        Endif      ! L_diff_ctl

        If ( CycleNo == NumCycles ) Then
          DEALLOCATE (BL_COEF_KM)
          DEALLOCATE (BL_COEF_KH)
          DEALLOCATE (visc_m)
          DEALLOCATE (visc_h)
          DEALLOCATE (visc_BL_m)
          DEALLOCATE (RNEUTML)
          DEALLOCATE (shear)
          If (L_Physics)then
            DEALLOCATE (FM_3D)
            DEALLOCATE (FH_3D)
          End If     ! L_physics
        End If    !  CycleNo == NumCycles

      End If       ! ErrorStatus  ==  0

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS3 Diffusion',6)
!----------------------------------------------------------------------

      If ( L_tracer .and. CycleNo == NumCycles ) Then

! Obtain increments from Atmos_Physics2
! DEPENDS ON: tr_reset
        call TR_Reset(                                                  &
                      super_array_size, super_tracer_phys2,             &
                      L_CO2_interactive, CO2,                           &
                      L_Murk_advect, murk,                              &
                      L_Soot, soot_new,                                 &
                              soot_agd,                                 &
                              soot_cld,                                 &
                      L_SULPC_SO2, SO2,                                 &
                                   SO4_aitken,                          &
                                   so4_accu,                            &
                                   so4_diss,                            &
                      L_sulpc_nh3, nh3,                                 &
                      L_sulpc_dms, dms,                                 &
                      L_dust, DUST_DIV1,                                &
                              DUST_DIV2,                                &
                              DUST_DIV3,                                &
                              DUST_DIV4,                                &
                              DUST_DIV5,                                &
                              DUST_DIV6,                                &
                      L_biomass, bmass_new,                             &
                                 bmass_agd,                             &
                                 bmass_cld,                             &
                      L_ocff, ocff_new,                                 &
                              ocff_agd,                                 &
                              ocff_cld,                                 &
                      L_USE_CARIOLLE, OZONE_TRACER,                     &
                      tracer_phys2, tracer, tracer_ukca,                &
                      row_length, rows,                                 &
                      model_levels, tr_levels, tr_vars, tr_ukca,        &
                      offx, offy                                        &
                                            )

      end If ! L_tracer .and. CycleNo == NumCycles

! Obtain increments from Atmos_Physics2 and diffusion

      If (L_moist_nonhydro_conserve)then
! after atmos_phys2 q_star holds q_dash (latest estimate to q_(n+1))
        IF(l_mix_ratio)then
          allocate ( mix_v_inter(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
          allocate ( mix_cl_inter(1-offx:row_length+offx,               &
     &                            1-offy:rows+offy, wet_levels) )
          allocate ( mix_cf_inter(1-offx:row_length+offx,               &
     &                            1-offy:rows+offy, wet_levels) )
          if(L_mcr_qcf2)then
            allocate ( mix_cf2_inter(1-offx:row_length+offx,            &
     &                               1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_cf2_inter(1,1,1) )
          endif
          if(L_mcr_qrain)then
            allocate ( mix_rain_inter(1-offx:row_length+offx,           &
     &                                1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_rain_inter(1,1,1) )
          endif
          if(L_mcr_qgraup)then
            allocate ( mix_graup_inter(1-offx:row_length+offx,          &
     &                                 1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_graup_inter(1,1,1) )
          endif
! DEPENDS ON: q_to_mix
          call q_to_mix (row_length, rows, wet_levels,                  &
     &                   offx,offy     ,                                &
     &                   q_star, qcl_star, qcf_star,                    &
     &               qcf2_star, qrain_star, qgraup_star,                &
!    &               .false.,.false.,.false.,
     &                L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,            &
     &               mix_v_inter, mix_cl_inter, mix_cf_inter,           &
     &               mix_cf2_inter, mix_rain_inter, mix_graup_inter     &
     &                   )
          do k=1, wet_levels
            do j=1,rows
              do i=1,row_length
                mix_v_phys2(i,j,k) = mix_v_inter(i,j,k)                 &
     &                             - mix_v_phys2(i,j,k)
                mix_cl_phys2(i,j,k) = mix_cl_inter(i,j,k)               &
     &                              - mix_cl_phys2(i,j,k)
                mix_cf_phys2(i,j,k) = mix_cf_inter(i,j,k)               &
     &                              - mix_cf_phys2(i,j,k)
              end do
            end do
          end do
          if(L_mcr_qcf2)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  mix_cf2_phys2(i,j,k) = mix_cf2_inter(i,j,k)           &
     &                                 - mix_cf2_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_mcr_qrain)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  mix_rain_phys2(i,j,k) = mix_rain_inter(i,j,k)         &
     &                                  - mix_rain_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_mcr_qgraup)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  mix_graup_phys2(i,j,k) = mix_graup_inter(i,j,k)       &
     &                                   - mix_graup_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_pc2)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  cf_phys2(i,j,k) = cf_star(i,j,k)                      &
     &                            - cf_phys2(i,j,k)
                  cfl_phys2(i,j,k) = cfl_star(i,j,k)                    &
     &                             - cfl_phys2(i,j,k)
                  cff_phys2(i,j,k) = cff_star(i,j,k)                    &
     &                             - cff_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          deallocate (mix_v_inter)
          deallocate (mix_cl_inter)
          deallocate (mix_cf_inter)
          deallocate (mix_cf2_inter)
          deallocate (mix_rain_inter)
          deallocate (mix_graup_inter)

        else

          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_phys2(i,j,k) = q_star(i,j,k)                          &
     &                         - q_phys2(i,j,k)
                qcl_phys2(i,j,k) = qcl_star(i,j,k)                      &
     &                           - qcl_phys2(i,j,k)
                qcf_phys2(i,j,k) = qcf_star(i,j,k)                      &
     &                           - qcf_phys2(i,j,k)
              End Do
            End Do
          End Do
          if(L_pc2)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  cf_phys2(i,j,k) = cf_star(i,j,k)                      &
     &                            - cf_phys2(i,j,k)
                  cfl_phys2(i,j,k) = cfl_star(i,j,k)                    &
     &                             - cfl_phys2(i,j,k)
                  cff_phys2(i,j,k) = cff_star(i,j,k)                    &
     &                             - cff_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_mcr_qcf2)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  qcf2_phys2(i,j,k) = qcf2_star(i,j,k)                  &
     &                              - qcf2_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_mcr_qrain)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  qrain_phys2(i,j,k) = qrain_star(i,j,k)                &
     &                               - qrain_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_mcr_qgraup)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  qgraup_phys2(i,j,k) = qgraup_star(i,j,k)              &
     &                               -  qgraup_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
        endif !L_mix_ratio

      end if  !L_moist_nonhydro_conserve

! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after diffusion  ***'
        ENDIF ! L_print_pe .or. mype ==0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      endif !  L_print_L2norms
! ----------------------------------------------------------------------
! Section 4.  Form and solve Helmholtz equation and update variables.
! ----------------------------------------------------------------------
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AS4 Solver',5)

      If (ErrorStatus  ==  0) then

! Apply diagnostics only at last cycle.
        If ( CycleNo == NumCycles ) Then

! Save current values to form diagnostics increments over Solver
        IF( sf(185,10)) THEN
          ALLOCATE ( u_incr_diagnostic(row_length,rows,model_levels) )
! Hold u increment (no halos needed)
          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                u_incr_diagnostic(i,j,k) =R_u(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                   ! on STASHflag

        IF( sf(186,10) ) THEN
          ALLOCATE ( v_incr_diagnostic(row_length,n_rows,model_levels) )
! Hold v increment (no halos needed)
          DO k=1,model_levels
            DO j=1,n_rows
              DO i=1,row_length
                v_incr_diagnostic(i,j,k) =R_v(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(187,10) ) THEN
          ALLOCATE ( w_incr_diagnostic(row_length,rows,model_levels) )
! Hold w increment (no halos needed)
          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                w_incr_diagnostic(i,j,k) =R_w(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        End If ! CycleNo == NumCycles
! set halos for q_star, qcl_star, qcf_star and theta_star

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   q_star, row_length, rows, wet_levels,          &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(q_star,row_length,rows,                &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   qcl_star, row_length, rows, wet_levels,        &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(qcl_star,row_length,rows,              &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   qcf_star, row_length, rows, wet_levels,        &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(qcf_star,row_length,rows,              &
     &                           wet_levels,Offx,Offy)
        If(L_mcr_qcf2)then
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
     &                   qcf2_star, row_length, rows, wet_levels,       &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(qcf2_star,row_length,rows,           &
     &                           wet_levels,Offx,Offy)
        endif
        If(L_mcr_qrain)then
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
     &                   qrain_star, row_length, rows, wet_levels,      &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(qrain_star,row_length,rows,          &
     &                           wet_levels,Offx,Offy)
        endif
        If(L_mcr_qgraup)then
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
     &                   qgraup_star, row_length, rows, wet_levels,     &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(qgraup_star,row_length,rows,         &
     &                           wet_levels,Offx,Offy)
        endif


! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   theta_star, row_length, rows, model_levels,    &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(theta_star,row_length,rows,            &
     &                           model_levels,Offx,Offy)

        If (L_mix_ratio) Then

! DEPENDS ON: q_to_mix_halo
        call q_to_mix_halo (row_length, rows, wet_levels,               &
     &                 offx,offy     ,                                  &
     &                 q_star, qcl_star, qcf_star,                      &
     &                 qcf2_star,qrain_star, qgraup_star,               &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 mix_v_star, mix_cl_star, mix_cf_star,            &
     &                 mix_cf2_star, mix_rain_star, mix_graup_star      &
     &                 )

        End If       !  L_mix_ratio

! Calculate updated values of pressure gradient terms.
! DEPENDS ON: timer
        If (Ltimer) Call timer ('pg_update',3)

! DEPENDS ON: ni_pg_update
        Call NI_pg_update(                                              &
     &                 Q, QCL, QCF,                                     &
     &                 QCF2, QRAIN, QGRAUP,                             &
     &                 q_star, qcl_star, qcf_star,                      &
     &                 qcf2_star, qrain_star, qgraup_star,              &
     &                 mix_v, mix_cl, mix_cf,                           &
     &                 mix_cf2, mix_rain, mix_graup,                    &
     &                 mix_v_star, mix_cl_star, mix_cf_star,            &
     &                 mix_cf2_star, mix_rain_star, mix_graup_star,     &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 theta_star, THETA, EXNER_RHO_LEVELS,             &
     &                 r_theta_levels, r_rho_levels,                    &
     &                 sec_theta_latitude,                              &
     &                 delta_lambda, delta_phi, timestep,               &
     &                 grecip_dlamp(lambda_start), recip_dphip,         &
     &                 wt_lambda_p, wt_phi_p,                           &
     &                 alpha_3, alpha_4,                                &
     &                 Cp, epsilon,                                     &
     &                 row_length, rows, n_rows, model_levels,          &
     &                 wet_levels, model_domain,                        &
     &                 first_constant_r_rho_level,                      &
     &                 offx, offy, halo_i, halo_j, at_extremity,        &
     &                 CycleNo, L_new_tdisc, L_qwaterload,              &
     &                 R_u, R_v, R_w, L_mix_ratio, L_regular )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('pg_update',4)

! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after pg_update  ***'
        ENDIF ! L_print_pe .or. mype ==0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      endif !  L_print_L2norms
! ---------------------------------------------------------------

        If (model_domain  ==  mt_global ) Then
          If ( L_polar_filter_incs .or. L_filter_incs) Then

! DEPENDS ON: ni_filter_incs_ctl
            Call NI_filter_incs_Ctl(                                    &
     &                      THETA, theta_star, R_u, R_v, R_w,           &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      r_theta_levels, r_rho_levels,               &
     &                      r_at_u, r_at_v,                             &
     &                      delta_lambda, delta_phi,                    &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      cos_v_latitude, sec_v_latitude,             &
     &                      polar_filter_north_lat_limit,               &
     &                      polar_filter_south_lat_limit,               &
     &                      polar_filter_coefficient,                   &
     &                      polar_filter_n_sweeps,                      &
     &                      polar_filter_step_per_sweep,                &
     &                      polar_filter_lat_limit,                     &
     &                      max_121_rows, u_sweeps, v_sweeps,           &
     &                      global_u_filter, global_v_filter,           &
     &                      u_begin, u_end, v_begin, v_end,             &
     &                      diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
     &                      diff_coeff_thermo, diff_coeff_wind,         &
     &                      diff_order_thermo, diff_order_wind,         &
     &                      horizontal_level, mype,                     &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      nproc, nproc_x, nproc_y, datastart,         &
     &                      neighbour, at_extremity, model_domain,      &
     &                      gc_proc_row_group, gc_proc_col_group,       &
     &                      L_polar_filter_incs, L_filter_incs,         &
     &                      L_pfcomb, L_pftheta, L_pfuv, L_pfw,         &
     &                      L_pofil_new, L_diff_incs, Ltimer,           &
     &                      exner_theta_levels,                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &                      STASHwork13)

            IF( SF(0,13) ) THEN  ! Diagnostics required for this section

! DEPENDS ON: timer
              IF (Ltimer) CALL TIMER('STASH',3)

! DEPENDS ON: stash
              CALL STASH(a_sm,a_im,13,STASHwork13,                      &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &        ErrorStatus,Cmessage)

! DEPENDS ON: timer
              IF (Ltimer) CALL TIMER('STASH',4)

              DEALLOCATE (STASHwork13)

            ENDIF !   SF(0,13)

          End If   ! L_polar_filter_incs .or. L_filter_incs

! calculate R_u at the poles.

! DEPENDS ON: polar_vector_wind_n
          Call Polar_vector_wind_n(                                     &
     &                       R_v,                                       &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, model_levels, mag_vector_np,       &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       offx, offy, global_row_length,             &
     &                       gc_proc_row_group, at_extremity)

          If (at_extremity(PSouth) ) Then
            Do k = 1,model_levels
              Do i = 1, row_length
                R_u(i,1,k) = - mag_vector_sp(k) *                       &
     &                         sin ( lambda_a(i) - dir_vector_sp(k) )
              End Do
            End Do
          End If    ! at_extremity(PSouth)
          If (at_extremity(PNorth) ) Then
            Do k = 1,model_levels
              Do i = 1, row_length
                R_u(i,rows,k) = mag_vector_np(k) *                      &
     &                         sin ( lambda_a(i) - dir_vector_np(k) )
              End Do
            End Do
          End If    ! at_extremity(PNorth)

        End If    !model_domain

! Do halo swops for * variables

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   R_u, row_length, rows, model_levels,           &
     &                   offx, offy, fld_type_u, .true.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(R_u,row_length,rows,                   &
     &                           wet_levels,Offx,Offy)

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   R_v, row_length, n_rows, model_levels,         &
     &                   offx, offy, fld_type_v, .true.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(R_v,row_length,n_rows,                 &
     &                           wet_levels,Offx,Offy)

      End If        !   ErrorStatus  ==  0

! ----------------------------------------------------------------------
! Section 4.1 Form and solve Helmholtz equation, returning all
!             corrections that are required.
! ----------------------------------------------------------------------
      If (ErrorStatus  ==  0) then

! in the idealised versions the boundary layer coefficients Km are
! replaced by the frictional timescale.

        If (model_domain ==  mt_lam) Then

          L_do_halos=.FALSE.
          L_do_boundaries=.TRUE.

          If (RIM_STEPSA  ==  0) Then
            increment_factor=0.0
          Else
            increment_factor=1.0/                                       &
     &      (RIM_STEPSA-MOD(Timestep_Number-1,RIM_STEPSA))
          End If

          lbc_size=LENRIMA(fld_type_p,halo_type_extended,               &
     &             rima_type_norm)

          If (L_lbc_new) Then
! Obtain Exner tendency to pass to the solver to apply lbc

! Apply additional balance to EXNER_LBC_TEND, RHO_LBC_TEND 
! and W_LBC_TEND only on timestep after up_bound has been called.

          If ( MOD(Timestep_Number-1,RIM_STEPSA) == 0) Then 
            If (L_LBC_balance) Then

! DEPENDS ON: balance_lbc_values
              CALL BALANCE_LBC_VALUES(                                  &
     &        EXNER_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,             &
     &        Q_LBC_TEND, W_LBC_TEND, W_ADV_LBC_TEND,                   &
     &        U_LBC_TEND, V_LBC_TEND,                                   &
     &        R_RHO_LEVELS, R_THETA_LEVELS,                             &
     &        ROW_LENGTH, ROWS, model_levels, wet_levels, HALO_I,HALO_J,&
     &        LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),    &
     &        LENRIMA(fld_type_u,halo_type_extended,rima_type_norm),    &
     &        LENRIMA(fld_type_v,halo_type_extended,rima_type_norm),    &
     &        LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),&
     &        LBC_STARTA(1,fld_type_u,halo_type_extended,rima_type_norm),&
     &        LBC_STARTA(1,fld_type_v,halo_type_extended,rima_type_norm),&
     &        RIMWIDTHA, N_RIMS_TO_DO, RIMWEIGHTSA, AT_EXTREMITY,       &
     &        DELTA_PHI, DELTA_LAMBDA,                                  &
     &        BASE_PHI, BASE_LAMBDA,                                    &
     &        TWO_OMEGA, DATASTART,                                     &
     &        LAT_ROT_NP,                                               &
     &        GLOBAL_ROW_LENGTH, GLOBAL_ROWS, L_int_uvw_lbc             &
     &        )

            End If  !  L_LBC_balance
          End If  !  MOD(Timestep_Number-1,RIM_STEPSA) == 0

          DO k = 1, MODEL_LEVELS
            DO i=1,lbc_size
              EXNER_LBC_REAL_TEND(i,k) = increment_factor *             &
     &                         ( EXNER_LBC_TEND(i,k) - EXNER_LBC(i,k) )
            END DO
          END DO

          Else  ! Original lbc lgorithm requires winds at boundaries

! If old lbcs set the outer boundaries of R_u, R_v and R_w to the 
! difference between the boundary at time level n+1 and time level n

          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              W_LBC_REAL_TEND(i,k) = increment_factor *                 &
     &             (W_LBC_TEND(i,k) - W_LBC(i,k))
            END DO
          END DO

          lbc_size=LENRIMA(fld_type_u,halo_type_extended,               &
     &             rima_type_norm)
          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              U_LBC_REAL_TEND(i,k) = increment_factor *                 &
     &             (U_LBC_TEND(i,k) - U_LBC(i,k))
            END DO
          END DO
          lbc_size=LENRIMA(fld_type_v,halo_type_extended,               &
     &             rima_type_norm)

          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              V_LBC_REAL_TEND(i,k) = increment_factor *                 &
     &             (V_LBC_TEND(i,k) - V_LBC(i,k))
            END DO
          END DO
! U
! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH,ROWS,Offx,Offy,                                    &
     &    MODEL_LEVELS,fld_type_u,R_u,                                  &
     &    LENRIMA(fld_type_u,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_u,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_u,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    U_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm),n_rims_to_do,                       &
     &    RIMWEIGHTSA,AT_EXTREMITY,                                     &
     &    L_do_boundaries,L_do_halos)

! V
! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH,N_ROWS,Offx,Offy,                                  &
     &    MODEL_LEVELS,fld_type_v,R_v,                                  &
     &    LENRIMA(fld_type_v,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_v,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_v,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    V_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm),n_rims_to_do,                       &
     &    RIMWEIGHTSA,AT_EXTREMITY,                                     &
     &    L_do_boundaries,L_do_halos)

! W
! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH,ROWS,0,0,                                          &
     &    MODEL_LEVELS,fld_type_p,R_w,                                  &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    W_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm),n_rims_to_do,                       &
     &    RIMWEIGHTSA,AT_EXTREMITY,                                     &
     &    L_do_boundaries,L_do_halos)

          End If !  L_lbc_new

        Else If (model_domain  ==  mt_cyclic_lam) Then
! store old time level u and v in R_u and R_v at boundaries
          If (at_extremity(PSouth) ) Then
            Do k = 1, model_levels
              Do i = 1-offx, row_length+offx
                R_u(i,1,k)=U(i,1,k)
                R_v(i,1,k)=V(i,1,k)
              End Do
            End Do
            Do k = 1, model_levels
              Do i = 1, row_length
                R_w(i,1,k)=W(i,1,k)
              End Do
            End Do
          End If    ! at_extremity(PSouth)
          If (at_extremity(PNorth) ) Then
            Do k = 1, model_levels
              Do i = 1-offx, row_length+offx
                R_u(i,rows,k)= U(i,rows,k)
                R_v(i,n_rows,k)=V(i,n_rows,k)
              End Do
            End Do
            Do k = 1, model_levels
              Do i = 1, row_length
                R_w(i,rows,k)=W(i,rows,k)
              End Do
            End Do
          End If    ! at_extremity(PNorth)

! set increments to u, v and w at boundaries.
          If (at_extremity(PSouth) ) Then

            Do k = 1, model_levels
              Do i = 1-offx, row_length+offx
                R_u(i,1,k)=U(i,1,k)-R_u(i,1,k)
                R_v(i,1,k)=V(i,1,k)-R_v(i,1,k)
              End Do
            End Do
            Do k = 1, model_levels
              Do i = 1, row_length
                R_w(i,1,k)=W(i,1,k)-R_w(i,1,k)
              End Do
            End Do
          End If    ! at_extremity(PSouth)
          If (at_extremity(PNorth) ) Then
            Do k = 1, model_levels
              Do i = 1-offx, row_length+offx
                R_u(i,rows,k) = U(i,rows,k) - R_u(i,rows,k)
                R_v(i,n_rows,k) = V(i,n_rows,k) - R_v(i,n_rows,k)
              End Do
            End Do
            Do k = 1, model_levels
              Do i = 1, row_length
                R_w(i,rows,k) = W(i,rows,k) - R_w(i,rows,k)
              End Do
            End Do
          End If    ! at_extremity(PNorth)

        Endif  ! model_domain

! Note:
! R_u, R_v, R_w on output contain u_prime, v_prime, w_prime

! DEPENDS ON: timer
        If (Ltimer) Call timer ('PE_Helmholtz',3)

        ALLOCATE ( dtheta_dr_term(row_length,rows,model_levels) )


        If ( CycleNo == 1 ) Then
          GCR_zero_guess_it = GCR_zero_init_guess
        Else
          GCR_zero_guess_it = ( .not. L_GCR_cycle_opt )                 &
     &                          .and. GCR_zero_init_guess
          If ( .NOT. L_GCR_cycle_opt ) Then
              exner_prime(:,:,:) = 0.0
          End If
        End If

        If ( CycleNo == 1 ) Then
          If ( GCR_use_tol_abs ) Then
            GCR_run_tol_abs = GCR_tol_abs
            GCR_run_tol_res = 0.0
          Else
            GCR_run_tol_res = GCR_tol_res
            GCR_run_tol_abs = 0.0
          End If
        Else
          If ( GCR_use_tol_abs ) Then
            GCR_run_tol_abs = GCR_tol_abs2
            GCR_run_tol_res = 0.0
          Else
            GCR_run_tol_res = GCR_tol_res2
            GCR_run_tol_abs = 0.0
          End If
        End If

! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms before solver  ***'
        ENDIF ! L_print_pe .or. mype ==0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      endif !  L_print_L2norms
! ---------------------------------------------------------------

! DEPENDS ON: ni_pe_helmholtz
        Call NI_PE_Helmholtz(                                           &
     &                      U, V, W,                                    &
     &                      r_theta_levels,                             &
     &                      r_rho_levels, P, RHO,                       &
     &                      rho_np1, THETA,                             &
     &                      theta_star,theta_np1,                       &
     &                      Q, QCL, QCF,                                &
     &                      QCF2, QRAIN, QGRAUP,                        &
     &                      mix_v, mix_cl, mix_cf,                      &
     &                      mix_cf2, mix_rain, mix_graup,               &
     &                      mix_v_star, mix_cl_star, mix_cf_star,       &
     &                      mix_cf2_star, mix_rain_star, mix_graup_star,&
     &                      mix_v_np1, mix_cl_np1, mix_cf_np1,          &
     &                      mix_cf2_np1, mix_rain_np1, mix_graup_np1,   &
     &                      q_star, q_np1, qcl_star, qcf_star,          &
     &                      qcf2_star, qrain_star, qgraup_star,         &
     &                      qcl_np1, qcf_np1,                           &
     &                      qcf2_np1 , qrain_np1, qgraup_np1,           &
     &                      rho_Km, cH, G_term_tol,                     &
     &                      EXNER_RHO_LEVELS,                           &
     &                      EXNER_THETA_LEVELS,                         &
     &                      frictional_timescale,                       &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      cos_v_latitude, sec_v_latitude,             &
     &                      FV_cos_theta_latitude,                      &
     &                      FV_sec_theta_latitude,                      &
     &                      f3_at_u, f3_at_v,                           &
     &                      timestep, timestep_number,                  &
     &                      row_length, rows, n_rows,                   &
     &                      model_levels, wet_levels,                   &
     &                      bl_levels, L_print_L2helm, L_print_pe,      &
     &                      L_flush6, norm_lev_start, norm_lev_end,     &
! ---------------------------------------------------------------
     &                      delta_lambda, delta_phi,                    &
     &                      glambda_p(lambda_start), phi_p,             &
     &                      glambda_u(lambda_start), phi_v,             &
     &                      gdlambda_p(lambda_start), dphi_p,           &
     &                      gdlambda_u(lambda_start), dphi_v,           &
     &                      grecip_dlamp(lambda_start), recip_dphip,    &
     &                      grecip_dlamu(lambda_start), recip_dphiv,    &
     &                      wt_lambda_p, wt_phi_p,                      &
     &                      wt_lambda_u, wt_phi_v,                      &
     &                      GCR_max_iterations, GCR_diagnostics,        &
     &                      GCR_its_switch(CycleNo), GCR_its_avg_step,  &
     &                      GCR_max_its(CycleNo), GCR_min_its(CycleNo), &
     &                      GCR_sum_its(CycleNo), GCR_max_time(CycleNo),&
     &                      GCR_min_time(CycleNo),                      &
     &  GCR_run_tol_res, GCR_run_tol_abs, GCR_use_tol_abs,              &
     &                      GCR_zero_guess_it,                          &
     &                      GCR_use_residual_Tol,                       &
     &                      GCR_adi_add_full_soln, L_gcr_fast_x,        &
     &                      GCR_precon_option, GCR_ADI_Pseudo_timestep, &
     &                      GCR_n_ADI_pseudo_timesteps,                 &
     &                      eta_theta_levels, eta_rho_levels,           &
     &                      alpha1, alpha2, alpha3, alpha4,             &
     &                      alpha_Cd, kappa, Cp, R, Pi,                 &
     &                      epsilon, model_domain, L_physics,           &
     &                      GCR_Restart_value,                          &
     &                      first_constant_r_rho_level,                 &
     &                      first_constant_r_rho_level_m1,              &
     &                      R_u, R_v, R_w, exner_prime, dtheta_dr_term, &
     &                      EXNER_LBC_REAL_TEND,                        &
     &                      LENRIMA(1,1,rima_type_norm),                &
     &                      LBC_SIZEA(1,1,1,rima_type_norm),            &
     &                      LBC_STARTA(1,1,1,rima_type_norm),           &
     &                      RIMWIDTHA(rima_type_norm), RIMWEIGHTSA,     &
!!! parallel variables
     &                      mype, nproc, nproc_x, nproc_y,              &
     &                      halo_i, halo_j, datastart,                  &
     &                      L_regular, at_extremity,                    &
     &                      n_rims_to_do, offx, offy,                   &
     &                      gc_proc_row_group, gc_proc_col_group,       &
     &                      global_row_length, global_rows,             &
     &                      g_rows, g_row_length, g_datastart, CycleNo, &
     &                      L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
     &                      L_mix_ratio, L_fint_theta, L_new_tdisc,     &
     &                      L_qwaterload, L_lbc_new                     &
     &                      )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('PE_Helmholtz',4)

! If increments are swapped now then this avoids needing swap
! in update rho

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   R_u, row_length, rows, model_levels,           &
     &                   offx, offy, fld_type_u, .true.)

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   R_v, row_length, n_rows, model_levels,         &
     &                   offx, offy, fld_type_v, .true.)


!----------------------------------------------------------------------
! Section 4.1 Diagnostics from Solver (may need to move to after
!  rho update if wanted to output rho increments in future).
!----------------------------------------------------------------------
! Apply diagnostics only at lat cycle.
        If ( CycleNo == NumCycles ) Then

! DEPENDS ON: timer
        if(Ltimer) Call timer ('Diag_sol',3)
!
        IF (SF(0,10)                                                    &
                         ! diagnostics required for this section
     &      .and. ErrorStatus == 0) THEN

! Allocate diagnostic space for STASH
          ALLOCATE (STASHwork10(STASH_maxlen(10,A_im)))

! DEPENDS ON: diagnostics_solver
          CALL Diagnostics_solver(                                      &
     &       row_length,rows,n_rows,model_levels,                       &
! wind field increments after  solver :
     &       R_u,R_v,R_w,                                               &
! wind field increments before solver (on stashflag):
     &       u_incr_diagnostic,v_incr_diagnostic,w_incr_diagnostic,     &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &       STASHwork10)

! Tidy allocatable arrays
          IF( sf(185,10) ) THEN
            DEALLOCATE ( u_incr_diagnostic )
          ENDIF  ! on STASHflag
          IF( sf(186,10) ) THEN
            DEALLOCATE ( v_incr_diagnostic )
          ENDIF  ! on STASHflag
          IF( sf(187,10) ) THEN
            DEALLOCATE ( w_incr_diagnostic )
          ENDIF  ! on STASHflag

! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)

! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,10,STASHwork10,                          &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)

          DEALLOCATE (STASHwork10)

! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)

        ENDIF !   SF(0,10)

! DEPENDS ON: timer
        if(Ltimer) Call timer ('Diag_sol',4)

        End If ! CycleNo == NumCycles

      End If  !      ErrorStatus  ==  0

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AS4 Solver',6)

! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after solver  ***'
        ENDIF ! L_print_pe .or. mype ==0
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
      endif !  L_print_L2norms
! ---------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 5.0 Calculate rho at new time level, using flux form of
!             equation.
! ----------------------------------------------------------------------
! Call timer for updates code
! DEPENDS ON: timer
      If (Ltimer .and. CycleNo == 1 ) CALL TIMER('AS5-8 Updates',5)

      If (ErrorStatus  ==  0) then

        If ( CycleNo == 1 ) Then
          If ( sf(188,30) ) Then
            ALLOCATE ( inc_rho( 1-offx:row_length+offx,                 &
     &                        1-offy:rows+offy, model_levels) )
          Else
            ALLOCATE ( inc_rho( 1,1,1) )
          End If
        End If

! stash diag 30188 only needed at last cycle.
        L_do_inc_rho = (CycleNo == NumCycles) .and. sf(188,30)

! Set rho_n at
        If ( CycleNo == NumCycles .and. ( L_tracer .or.                 &
     &       L_moist_nonhydro_conserve .or. L_do_inc_rho ) ) Then

! store rho_n before updating by flux_rho
! DEPENDS ON: copy_field
      CALL COPY_FIELD(RHO,rho_n                                         &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels, model_levels, 1, model_levels       &
     &,               Offx, Offy, Offx, Offy                            &
     &,               fld_type_p, .true., .false., .false.)

! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(rho_n,row_length, rows,              &
     &                     model_levels,offx,offy)

        End If ! ! CycleNo == NumCycles .and. ...


        If (  CycleNo == NumCycles .OR. L_new_tdisc ) Then

! DEPENDS ON: ni_update_rho
          Call NI_Update_Rho(                                           &
     &                 rho, rho_n, rho_np1, inc_rho,                    &
     &                 u, v, w, R_u, R_v, R_w,                          &
     &                 q, qcl,qcf,                                      &
     &                 qcf2, qrain, qgraup,                             &
     &                 q_star, qcl_star, qcf_star,                      &
     &                 qcf2_star, qrain_star, qgraup_star,              &
     &                 q_np1, qcl_np1, qcf_np1,                         &
     &                 qcf2_np1, qrain_np1, qgraup_np1,                 &
     &                 mix_v, mix_cl, mix_cf,                           &
     &                 mix_cf2, mix_rain, mix_graup,                    &
     &                 mix_v_star, mix_cl_star, mix_cf_star,            &
     &                 mix_cf2_star, mix_rain_star, mix_graup_star,     &
     &                 mix_v_np1, mix_cl_np1, mix_cf_np1,               &
     &                 mix_cf2_np1, mix_rain_np1, mix_graup_np1,        &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 timestep, CycleNo, NumCycles,                    &
     &                 rows, n_rows, row_length,                        &
     &                 model_levels, wet_levels, model_domain,          &
     &                 first_constant_r_rho_level,                      &
     &                 alpha_1, alpha_2, n_rims_to_do,                  &
     &                 nproc, gc_proc_row_group,                        &
     &                 L_regular, at_extremity, global_row_length,      &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 cos_v_latitude, delta_lambda, delta_phi,         &
     &                 gdlambda_p(lambda_start), dphi_p,                &
     &                 grecip_dlamu(lambda_start), recip_dphiv,         &
     &                 wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v,    &
     &                 r_theta_levels, r_rho_levels, eta_theta_levels,  &
     &                 eta_rho_levels, FV_sec_theta_latitude,           &
     &                 wet_to_dry_n, wet_to_dry_np1,                    &
     &                 L_do_inc_rho, L_new_tdisc,                       &
     &                 L_mix_ratio, L_dry )

        End If !  CycleNo == NumCycles .OR. L_new_tdisc

        If ( CycleNo == NumCycles ) Then

          if(L_mix_ratio .and. .not. L_moist_nonhydro_conserve)then
            deallocate(mix_v)
            deallocate(mix_cl)
            deallocate(mix_cf)
            deallocate(mix_cf2)
            deallocate(mix_rain)
            deallocate(mix_graup)
          endif

          If (NumCycles > 1 .and. L_mix_ratio                           &
     &                      .and. .not. L_moist_nonhydro_conserve) Then
            deallocate( mix_v_star_phys1)
            deallocate( mix_cl_star_phys1)
            deallocate( mix_cf_star_phys1)
            If ( L_mcr_qcf2 ) deallocate( mix_cf2_star_phys1)
            If ( L_mcr_qrain ) deallocate( mix_rain_star_phys1)
            If ( L_mcr_qgraup ) deallocate( mix_graup_star_phys1)
          End If

        If (L_tracer) then
! DEPENDS ON: timer
          If (Ltimer) CALL TIMER('SL_tracer2',3)

! DEPENDS ON: sl_tracer2
          Call SL_tracer2(                                              &
     &                           super_array_size,                      &
     &                  super_tracer_phys1,super_tracer_phys2,          &
     &                 eta_theta_levels,                                &
     &                 r_rho_levels,  r_theta_levels,                   &
     &                 rho_n, RHO,                                      &
     &                 row_length, rows, model_levels,                  &
     &                 delta_lambda, delta_phi,                         &
     &                 glambda_p, phi_p, grecip_dlamp, recip_dphip,     & 
     &                 lambda_p_rm, lambda_p_rp, phi_p_rm, phi_p_rp,    &
     &                 recip_lambda_p_m, recip_lambda_p_0,              &
     &                 recip_lambda_p_p, recip_lambda_p_p2,             &
     &                 recip_phi_p_m, recip_phi_p_0,                    &
     &                 recip_phi_p_p, recip_phi_p_p2,                   &
     &                 base_lambda, base_phi,                           &
     &                 recip_dlam, recip_dphi, max_look,                &
     &                 look_lam, look_phi, halo_lam, halo_phi,          &
     &                 FV_cos_theta_latitude,                           &
     &                 wet_to_dry_n, wet_to_dry_np1,                    &
     &                 mype, nproc, nproc_x, nproc_y,                   &
     &                 halo_i, halo_j, datastart,                       &
     &                 g_i_pe, at_extremity,                            &
     &                 global_row_length,                               &
     &                 global_rows,                                     &
     &                 gc_all_proc_group,                               &
     &                 gc_proc_row_group,                               &
     &                 gc_proc_col_group, offx, offy,                   &
     &                 L_regular, L_sl_halo_reprod,                     &
     &                 high_order_scheme(moist_SL),                     &
     &                 monotone_scheme(moist_SL),                       &
     &                 model_domain, L_high(moist_SL),                  &
     &                 L_mono(moist_SL),                                &
     &                 L_conserve_tracers,                              &
     &                 check_bottom_levels,                             &
     &                 interp_vertical_search_tol,                      &
     &                 first_constant_r_rho_level,                      &
     &                 depart_lambda, depart_phi, depart_r_theta,       &
     &                 co2, L_CO2_interactive,                          &
     &                 Murk, L_Murk_advect,                             &
     &                 soot_new, soot_agd,                              &
     &                 soot_cld, L_soot,                                &
     &                 bmass_new, bmass_agd,                            &
     &                 bmass_cld, L_biomass,                            &
     &                 ocff_new, ocff_agd, ocff_cld, L_ocff,            &
     &                 DUST_DIV1,DUST_DIV2,                             &
     &                 DUST_DIV3,DUST_DIV4,                             &
     &                 DUST_DIV5,DUST_DIV6, L_DUST,                     &
     &                 so2, so4_aitken,                                 &
     &                 so4_accu,                                        &
     &                 so4_diss, nh3, dms,                              &
     &                 L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,           &
     &                 tracer, tr_levels, tr_vars,                      &
     &                 tracer_ukca, tr_ukca,                            &
     &                 tracer_phys1,tracer_phys2,                       &
     &                 i_start, i_stop, j_start, j_stop,                &
     &                 L_USE_CARIOLLE, OZONE_TRACER,                    &
     &                 ErrorStatus)

!      Cariolle scheme is called to calculate the tracer ozone. All the tracers
!      are calculated at the end of the timestep. The ozone tracer 
!      calculated here will be used in the radiation scheme on the next timestep.

!      Insert if statement here to check that tracer that goes into this call 
!      ozone. Ozone may not be the only tracer.
            If (PrintStatus  >=  PrStatus_Normal .AND.                  &
                first_atmstep_call) then

               WRITE(6,*) 'Atm_Step: L_USE_CARIOLLE = ',L_USE_CARIOLLE
            End If

            If (L_USE_CARIOLLE) then
               
               If (PrintStatus  >=  PrStatus_Normal .AND.                  &
                first_atmstep_call) then

                    WRITE(6,*) 'Atm_Step: Calling Cariolle_o3_psc'
               End If

! DEPENDS ON: cariolle_o3_psc
               CALL cariolle_o3_psc (OZONE_TRACER,                     &
                      O3_PROD_LOSS,   O3_P_L_VMR,                      &
                      O3_VMR,         O3_P_L_TEMP,                     &
                      O3_TEMP,        O3_P_L_COLO3,                    &
                      O3_COLO3,                                        &
                      THETA,                                           &
                      P_THETA_LEVELS,                                  &
                      offx,offy,theta_off_size,                        &
                      rows,row_length,                                 &
                      timestep,                                        &
                      EXNER_THETA_LEVELS,                              &
                      model_levels)      
!    
! Halos updated
! DEPENDS ON: swap_bounds
               call Swap_Bounds(OZONE_TRACER,                           &
                        row_length, rows, model_levels,                 &
                        offx, offy, fld_type_p, .false.)
!DEPENDS ON: fill_external_halos
               CALL FILL_EXTERNAL_HALOS(OZONE_TRACER,row_length, rows,  &
                         model_levels,offx,offy)
            
            else
               If (PrintStatus  >=  PrStatus_Normal .AND.               &
                   first_atmstep_call) then

                   WRITE(6,*) 'Atm_Step: Cariolle scheme not called'
               End If
            End if

! DEPENDS ON: timer
          If (Ltimer) CALL TIMER('SL_tracer2',4)

          If (.not. L_moist_nonhydro_conserve)then
            DEALLOCATE (depart_r_theta)
            DEALLOCATE (depart_lambda)
            DEALLOCATE (depart_phi)
          endif
        end if  ! L_tracer

! ----------------------------------------------------------------------
! Section 5.2 Recalculate sl_thermo to allow rho at n+1
!             to be used in the conservation step
!            Estimated values are returned in the _phys1 variables.
! ----------------------------------------------------------------------
!  At this point, q_star, qcl_star and qcf_star in the boundary layer
!  have been changed by the implicit bl calculation but are not used
!  from now on.
        If (L_moist_nonhydro_conserve) then

! DEPENDS ON: ni_sl_moist
          Call NI_SL_moist(                                             &
     &                 moisture_array_size,                             &
     &                 q, qcl, qcf,                                     &
     &                 cf_bulk,cf_frozen,                               &
     &                 cf_liquid,                                       &
     &                 qcf2, qrain, qgraup,                             &
     &                 q_phys1, qcl_phys1, qcf_phys1,                   &
     &                 cf_phys1, cff_phys1, cfl_phys1,                  &
     &                 qcf2_phys1, qrain_phys1, qgraup_phys1,           &
     &                 mix_v, mix_cl, mix_cf,                           &
     &                 mix_cf2, mix_rain, mix_graup,                    &
     &                 mix_v_phys1, mix_cl_phys1, mix_cf_phys1,         &
     &                 mix_cf2_phys1, mix_rain_phys1, mix_graup_phys1,  &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 L_pc2, L_mix_ratio,                              &
     &                 eta_theta_levels,                                &
     &                 r_rho_levels,  r_theta_levels,                   &
     &                 exner_theta_levels,                              &
     &                 rho_n, rho,                                      &
     &                 depart_lambda, depart_phi, depart_r_theta,       &
     &                 row_length, rows, n_rows,                        &
     &                 model_levels, wet_levels,                        &
     &                 delta_lambda, delta_phi,                         &
     &                 base_lambda, base_phi,                           &
     &                 glambda_p, phi_p, grecip_dlamp, recip_dphip,     &
     &                 lambda_p_rm, lambda_p_rp, phi_p_rm, phi_p_rp,    &
     &                 recip_lambda_p_m, recip_lambda_p_0,              &
     &                 recip_lambda_p_p, recip_lambda_p_p2,             &
     &                 recip_phi_p_m, recip_phi_p_0,                    &
     &                 recip_phi_p_p, recip_phi_p_p2,                   &
     &                 recip_dlam, recip_dphi, max_look,                &
     &                 look_lam, look_phi, halo_lam, halo_phi,          &
     &                 FV_cos_theta_latitude,                           &
     &                 wet_to_dry_n, wet_to_dry_np1,                    &
     &                 L_regular, n_rims_to_do,                         &
     &                 mype, nproc, nproc_x, nproc_y,                   &
     &                 halo_i, halo_j, datastart,                       &
     &                 g_i_pe, at_extremity,                            &
     &                 global_row_length, global_rows,                  &
     &                 gc_proc_row_group,                               &
     &                 gc_proc_col_group, offx, offy,                   &
     &                 L_sl_halo_reprod,                                &
     &                 high_order_scheme(moist_SL),                     &
     &                 monotone_scheme(moist_SL),                       &
     &                 model_domain, L_high(moist_SL),                  &
     &                 L_mono(moist_SL), L_conserv(moist_SL),           &
     &                 check_bottom_levels, interp_vertical_search_tol, &
     &                 first_constant_r_rho_level,                      &
     &                 ErrorStatus )

          DEALLOCATE (depart_lambda)
          DEALLOCATE (depart_phi)
          DEALLOCATE (depart_r_theta)
          if(L_mix_ratio)then
            deallocate(mix_v)
            deallocate(mix_cl)
            deallocate(mix_cf)
            deallocate(mix_cf2)
            deallocate(mix_rain)
            deallocate(mix_graup)
          endif

          If ( NumCycles > 1 .and. L_mix_ratio ) Then
            deallocate( mix_v_star_phys1)
            deallocate( mix_cl_star_phys1)
            deallocate( mix_cf_star_phys1)
            If ( L_mcr_qcf2 ) deallocate( mix_cf2_star_phys1)
            If ( L_mcr_qrain ) deallocate( mix_rain_star_phys1)
            If ( L_mcr_qgraup ) deallocate( mix_graup_star_phys1)
          End If
        end if          ! L_moist_nonhydro_conserve

      End If ! NoCycles == NumCycles

      End If       ! L_Primitive .and. ErrorStatus  ==  0

! ----------------------------------------------------------------------
! Section 6.0 Calculate u, v, w, theta, q, qcl, qcf at new time level.
!             Extrapolate u, v, w, to time level n+1.5 for use in
!             advection step on next timestep.
! ----------------------------------------------------------------------
      If ( ErrorStatus  ==  0) Then

        If ( CycleNo == 1 ) Then

        if(sf(0,30)) then
          if (sf(182,30) .or. sf(183,30) .or. sf(184,30)                &
     &       .or. sf(189,30) .or. sf(190,30) .or. sf(191,30)) then
            ALLOCATE ( inc_q(1-halo_i:row_length+halo_i,                &
     &                       1-halo_j:rows+halo_j,wet_levels) )
            ALLOCATE ( inc_qcl(1-halo_i:row_length+halo_i,              &
     &                         1-halo_j:rows+halo_j,wet_levels) )
            ALLOCATE ( inc_qcf(1-halo_i:row_length+halo_i,              &
     &                         1-halo_j:rows+halo_j,wet_levels) )
            ALLOCATE ( inc_cf(1-halo_i:row_length+halo_i,               &
     &                       1-halo_j:rows+halo_j,wet_levels) )
            ALLOCATE ( inc_cfl(1-halo_i:row_length+halo_i,              &
     &                         1-halo_j:rows+halo_j,wet_levels) )
            ALLOCATE ( inc_cff(1-halo_i:row_length+halo_i,              &
     &                         1-halo_j:rows+halo_j,wet_levels) )

!Store values for increment calculations for moisture variables
            Do k = 1, wet_levels
              Do j = j_start, j_stop
                Do i = i_start, i_stop
                  inc_q(i,j,k)   = q(i,j,k)
                  inc_qcl(i,j,k) = qcl(i,j,k)
                  inc_qcf(i,j,k) = qcf(i,j,k)
                  inc_cf(i,j,k)  = cf_bulk(i,j,k)
                  inc_cfl(i,j,k) = cf_liquid(i,j,k)
                  inc_cff(i,j,k) = cf_frozen(i,j,k)
                End Do
              End Do
            End Do
          else
            ALLOCATE ( inc_q(1,1,1) )
            ALLOCATE ( inc_qcl(1,1,1) )
            ALLOCATE ( inc_qcf(1,1,1) )
            ALLOCATE ( inc_cf(1,1,1) )
            ALLOCATE ( inc_cfl(1,1,1) )
            ALLOCATE ( inc_cff(1,1,1) )
          endif
        endif       ! sf(0,30)

        If ( sf(181,30) ) then
          ALLOCATE ( inc_t( 1-offx:row_length+offx,                     &
     &                      1-offy:rows+offy, model_levels) )
        else
          ALLOCATE ( inc_t( 1,1,1) )
        endif

        L_do_inc_vels = sf(185,30) .or. sf(186,30) .or. sf(187,30)

        if ( L_do_inc_vels ) then
          ALLOCATE ( inc_u( 1-offx:row_length+offx,                     &
     &                      1-offy:rows+offy, model_levels) )
          ALLOCATE ( inc_v( 1-offx:row_length+offx,                     &
     &                      1-offy:n_rows+offy, model_levels) )
          ALLOCATE ( inc_w(row_length, rows, model_levels) )
        else
          ALLOCATE ( inc_u( 1,1,1) )
          ALLOCATE ( inc_v( 1,1,1) )
          ALLOCATE ( inc_w( 1,1,1) )
        endif

        End If ! CycleNo == 1

! DEPENDS ON: update_fields
        call update_fields(                                             &
     &                   NumCycles, CycleNo, L_new_tdisc,               &
     &                   L_mix_ratio, extrp_weight,                     &
     &                   EXNER_RHO_LEVELS,                              &
     &                   EXNER_THETA_LEVELS,                            &
     &                   U, u_np1, V, v_np1,                            &
     &                   W, w_np1,                                      &
     &                   U_ADV, V_ADV, W_ADV,                           &
     &                   THETA,Q,QCL,QCF,                               &
     &                   QCF2, QRAIN, QGRAUP,                           &
     &                   CF_BULK,CF_LIQUID,                             &
     &                   CF_FROZEN,                                     &
     &                   exner_prime, R_u, R_v, R_w,                    &
     &                   theta_star, theta_np1, dtheta_dr_term,         &
     &                   q_star, q_np1, qcl_star, qcl_np1, qcf_star,    &
     &                   qcf_np1, qcf2_star, qcf2_np1, qrain_star,      &
     &                   qrain_np1, qgraup_star, qgraup_np1,            &
     &                   cf_star, cfl_star, cff_star,                   &
     &                   row_length, rows, n_rows, model_levels,        &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   sin_theta_longitude, cos_theta_longitude,      &
     &                   mag_vector_np, dir_vector_np,                  &
     &                   mag_vector_sp, dir_vector_sp,                  &
     &                   global_row_length, gc_proc_row_group,          &
     &                   at_extremity, datastart,model_domain,          &
     &                   alpha_2 , timestep, wet_levels,                &
     &                   i_start, i_stop, j_start, j_stop,              &
     &                   j_begin, j_end,                                &
     &                   delta_lambda, glambda_p(lambda_start),         &
     &                   glambda_u(lambda_start),                       &
     &                   inc_t, inc_u, inc_v, inc_w,                    &
     &                   sf(181,30), L_do_inc_vels, L_pc2,              &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   L_regular )

        If ( NumCycles > 1 .and. CycleNo < NumCycles ) Then
!
! For consistency in the formulation, the _adv and _np1 vars need
! to be updated in the lateral boundary region.
!
! DEPENDS ON: swap_bounds
          call Swap_Bounds(U_ADV,                                       &
     &                   row_length, rows, model_levels,                &
     &                   halo_i, halo_j, fld_type_u, .true.)

! DEPENDS ON: swap_bounds
          call Swap_Bounds(V_ADV,                                       &
     &                   row_length, n_rows, model_levels,              &
     &                   halo_i, halo_j, fld_type_v, .true.)

! DEPENDS ON: swap_bounds
          call Swap_Bounds(W_ADV,                                       &
     &                   row_length, rows, model_levels+1,              &
     &                   halo_i, halo_j, fld_type_p, .false.)

          If ( L_new_tdisc ) Then

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(u_np1,                                     &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, fld_type_u, .true.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(v_np1,                                     &
     &                   row_length, n_rows, model_levels,              &
     &                   offx, offy, fld_type_v, .true.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(w_np1,                                     &
     &                   row_length, rows, model_levels+1,              &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(                                           &
     &                   theta_np1, row_length, rows, model_levels,     &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(                                           &
     &                   rho_np1, row_length, rows, model_levels,       &
     &                   offx, offy, fld_type_p, .false.)


!
! NOTE: mix_v_star etc will temporarily hold q_star etc. These will be
!       converted to mix_v_np1 at the end of the current iteration.
!
            If ( L_mix_ratio ) Then

              mix_v_star = q_star
              mix_cl_star = qcl_star
              mix_cf_star = qcf_star
              If ( L_mcr_qcf2 ) mix_cf2_star = qcf2_star
              If ( L_mcr_qrain ) mix_rain_star = qrain_star
              If ( L_mcr_qgraup ) mix_graup_star = qgraup_star

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   mix_v_star, row_length, rows, wet_levels,      &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   mix_cl_star, row_length, rows, wet_levels,     &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   mix_cf_star, row_length, rows, wet_levels,     &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qcf2 )                                         &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   mix_cf2_star, row_length, rows, wet_levels,    &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qrain )                                        &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   mix_rain_star, row_length, rows, wet_levels,   &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qgraup )                                       &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   mix_graup_star, row_length, rows, wet_levels,  &
     &                   offx, offy, fld_type_p, .false.)

            Else

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   q_np1, row_length, rows, wet_levels,           &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   qcl_np1, row_length, rows, wet_levels,         &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   qcf_np1, row_length, rows, wet_levels,         &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qcf2 )                                         &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   qcf2_np1, row_length, rows, wet_levels,        &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qrain )                                        &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   qrain_np1, row_length, rows, wet_levels,       &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qgraup )                                       &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   qgraup_np1, row_length, rows, wet_levels,      &
     &                   offx, offy, fld_type_p, .false.)

            End If ! .NOT. L_mix_ratio

          End If ! If L_new_tdisc

! Now need to fix Uadv and temp n+1 dynamical vars at LB region.
          If ( model_domain == mt_lam                                   &
     &         .or.  model_domain == mt_cyclic_lam ) Then

            If ( L_new_tdisc ) Then
!
! Fill external haloes of theta_np1, q_np1,... following approach
! for theta_star, q_star.
!
               If ( .NOT. L_mix_ratio ) Then

! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(q_np1,row_length,rows,        &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(qcl_np1,row_length,rows,      &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(qcf_np1,row_length,rows,      &
     &                           wet_levels,Offx,Offy)

                 If ( L_mcr_qcf2 ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(qcf2_np1,row_length,rows,   &
     &                           wet_levels,Offx,Offy)
                 End If

                 If ( L_mcr_qrain ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(qrain_np1,row_length,rows,  &
     &                           wet_levels,Offx,Offy)
                 End If

                 If ( L_mcr_qgraup ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(qgraup_np1,row_length,rows, &
     &                           wet_levels,Offx,Offy)
                 End If

               Else
! mix_v_star etc temporarily hold q_star

! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(mix_v_star,row_length,rows,   &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(mix_cl_star,row_length,rows,  &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(mix_cf_star,row_length,rows,  &
     &                           wet_levels,Offx,Offy)
                 If ( L_mcr_qcf2 ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(mix_cf2_star,row_length,    &
     &                           rows,wet_levels,Offx,Offy)
                 End If
                 If ( L_mcr_qrain ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(mix_rain_star,row_length,   &
     &                           rows,wet_levels,Offx,Offy)
                 End If
                 If ( L_mcr_qgraup ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(mix_graup_star,row_length,  &
     &                           rows,wet_levels,Offx,Offy)
                 End If
               End If ! .NOT. L_mix_ratio

! DEPENDS ON: FILL_EXTERNAL_HALOS
               CALL FILL_EXTERNAL_HALOS(theta_np1,row_length,rows,      &
     &                         model_levels,Offx,Offy)

               If (RIM_STEPSA == 0) Then
                 increment_factor=0.0
               Else
                 increment_factor=1.0/                                  &
     &           (RIM_STEPSA-MOD(Timestep_Number-1,RIM_STEPSA))
               End If

               lbc_size=LENRIMA(fld_type_p,halo_type_extended,          &
     &                         rima_type_norm)
               DO k=1,MODEL_LEVELS
                 DO i=1,lbc_size
                   W_LBC_REAL_TEND(i,k)=W_LBC(i,k)                      &
     &          +  increment_factor*( W_LBC_TEND(i,k) - W_LBC(i,k) )
                 END DO
               END DO

               lbc_size=LENRIMA(fld_type_u,halo_type_extended,          &
     &                          rima_type_norm)
               DO k=1,MODEL_LEVELS
                 DO i=1,lbc_size
                   U_LBC_REAL_TEND(i,k)=U_LBC(i,k)                      &
     &        + increment_factor*( U_LBC_TEND(i,k) - U_LBC(i,k) )      
                 END DO
               END DO

               lbc_size=LENRIMA(fld_type_v,halo_type_extended,          &
     &                         rima_type_norm)
               DO k=1,MODEL_LEVELS
                 DO i=1,lbc_size
                   V_LBC_REAL_TEND(i,k)=V_LBC(i,k)                      &
     &        + increment_factor*( V_LBC_TEND(i,k) - V_LBC(i,k) )
                 END DO
               END DO

               If ( CycleNo == 1 ) Then
                 Allocate (                                             &
     &                    RHO_LBC_REAL_TEND(LENRIMA(fld_type_p,         &
     &                    halo_type_extended,rima_type_norm),           &
     &                    MODEL_LEVELS) )
               End If

               lbc_size=LENRIMA(fld_type_p,halo_type_extended,          &
     &                         rima_type_norm)
               DO k=1,MODEL_LEVELS
                 DO i=1,lbc_size
                 RHO_LBC_REAL_TEND(i,k)=RHO_LBC(i,k)                    &
     &        + increment_factor*(RHO_LBC_TEND(i,k) - RHO_LBC(i,k) )
                 END DO
               END DO

            Else

              If ( CycleNo == 1 ) Allocate ( RHO_LBC_REAL_TEND(1,1) )

            End If ! If L_new_tdisc

! DEPENDS ON: TIMER
            IF (LTIMER) CALL TIMER('UPDATE_LBC_ITERSL',3)

            If ( .NOT. L_mix_ratio ) Then
! DEPENDS ON: UPDATE_LBC_ITERSL
              CALL UPDATE_LBC_ITERSL(                                   &
     &           ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,        &
     &           OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY, L_new_tdisc,     &
     &           L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                       &
     &           L_mcr_qgraup_lbc, L_pc2_lbc,                           &
     &           RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                 &
     &           LENRIMA(1,1,rima_type_norm),                           &
     &           LBC_SIZEA(1,1,1,rima_type_norm),                       &
     &           LBC_STARTA(1,1,1,rima_type_norm),                      &
     &           THETA_LBC,Q_LBC,QCL_LBC,                               &
     &           QCF_LBC,QCF2_LBC,QRAIN_LBC,                            &
     &           QGRAUP_LBC, CF_BULK_LBC,CF_LIQUID_LBC,                 &
     &           CF_FROZEN_LBC, RHO_LBC_REAL_TEND, EXNER_LBC,           &
     &           U_LBC_REAL_TEND,V_LBC_REAL_TEND,W_LBC_REAL_TEND,       &
     &           U_ADV_LBC,V_ADV_LBC,W_ADV_LBC,                         &
     &           theta_np1, q_np1, qcl_np1, qcf_np1,                    &
     &           qcf2_np1, qrain_np1, qgraup_np1,                       &
     &           CF_BULK,CF_LIQUID,CF_FROZEN,                           &
     &           rho_np1,EXNER_RHO_LEVELS,                              &
     &           u_np1,v_np1,w_np1,                                     &
     &           U_ADV,V_ADV,W_ADV)          
            Else
! mix_v_star etc hold q_star etc temporarily
! DEPENDS ON: UPDATE_LBC_ITERSL
              CALL UPDATE_LBC_ITERSL(                                   &
     &           ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,        &
     &           OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY, L_new_tdisc,     &
     &           L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                       &
     &           L_mcr_qgraup_lbc, L_pc2_lbc,                           &
     &           RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                 &
     &           LENRIMA(1,1,rima_type_norm),                           &
     &           LBC_SIZEA(1,1,1,rima_type_norm),                       &
     &           LBC_STARTA(1,1,1,rima_type_norm),                      &
     &           THETA_LBC,Q_LBC,QCL_LBC,                               &
     &           QCF_LBC,QCF2_LBC,QRAIN_LBC,                            &
     &           QGRAUP_LBC, CF_BULK_LBC,CF_LIQUID_LBC,                 &
     &           CF_FROZEN_LBC, RHO_LBC_REAL_TEND, EXNER_LBC,           &
     &           U_LBC_REAL_TEND,V_LBC_REAL_TEND,W_LBC_REAL_TEND,       &
     &           U_ADV_LBC,V_ADV_LBC,W_ADV_LBC,                         &
     &           theta_np1, mix_v_star, mix_cl_star, mix_cf_star,       &
     &           mix_cf2_star, mix_rain_star, mix_graup_star,           &
     &           CF_BULK,CF_LIQUID,CF_FROZEN,                           &
     &           rho_np1,EXNER_RHO_LEVELS,                              &
     &           u_np1,v_np1,w_np1,                                     &
     &           U_ADV,V_ADV,W_ADV)
            End If ! .NOT. L_mix_ratio

! DEPENDS ON: TIMER
            IF (LTIMER) CALL TIMER('UPDATE_LBC_ITERSL',4)

            If ( CycleNo == NumCycles-1 ) Deallocate(RHO_LBC_REAL_TEND)

          End If ! If model_domain == mt_lam ...

! Now convert qX(=mix_v_star,...) to mX and store at mix_v_np1 etc
          If ( L_mix_ratio .and. L_new_tdisc ) Then
            call q_to_mix_halo (row_length, rows, wet_levels,           &
     &                 offx, offy, mix_v_star, mix_cl_star, mix_cf_star,&
     &                 mix_cf2_star, mix_rain_star, mix_graup_star,     &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 mix_v_np1, mix_cl_np1, mix_cf_np1,               &
     &                 mix_cf2_np1, mix_rain_np1, mix_graup_np1         &
     &                 )
          End If

        End If ! NumCycles > 1 .and. CycleNo < NumCycles
!  LBC updating now done at end of time step     
! ---------------------------------------------------------------------
!   Section 6.1  Update lbcs for LAMs
! ---------------------------------------------------------------------

      If ( model_domain == mt_lam ) Then
      
        IF (L_Fixed_lbcs) Then
          L_update_lbcs = .false.
          L_apply_lbcs = .true.
        Else IF ( RIM_STEPSA == 0 ) Then
          L_update_lbcs = .false.
        Else IF ( L_lbc_new ) Then
          L_update_lbcs = .true.
        Else
          L_update_lbcs = .false.
        END IF !  L_Fixed_lbcs

        If ( L_update_lbcs ) Then
        
! DEPENDS ON: BOUNDVAL
              Call BOUNDVAL(LENRIMA(1,1,rima_type_norm),                &
     &                      L_mcr_qcf2_lbc, L_mcr_qrain_lbc,            &
     &                      L_mcr_qgraup_lbc, L_pc2_lbc,                &
     &                      L_murk_lbc, L_int_uvw_lbc,                  &
     &                      U_LBC, U_LBC_TEND,                          &
     &                      V_LBC, V_LBC_TEND,                          &
     &                      W_LBC, W_LBC_TEND,                          &
     &                      RHO_LBC, RHO_LBC_TEND,                      &
     &                      THETA_LBC, THETA_LBC_TEND,                  &
     &                      Q_LBC, Q_LBC_TEND,                          &
     &                      QCL_LBC, QCL_LBC_TEND,                      &
     &                      QCF_LBC, QCF_LBC_TEND,                      &
     &                      QCF2_LBC, QCF2_LBC_TEND,                    &
     &                      QRAIN_LBC, QRAIN_LBC_TEND,                  &
     &                      QGRAUP_LBC, QGRAUP_LBC_TEND,                &
     &                      CF_BULK_LBC, CF_BULK_LBC_TEND,              &
     &                      CF_LIQUID_LBC, CF_LIQUID_LBC_TEND,          &
     &                      CF_FROZEN_LBC, CF_FROZEN_LBC_TEND,          &
     &                      EXNER_LBC, EXNER_LBC_TEND,                  &
     &                      U_ADV_LBC, U_ADV_LBC_TEND,                  &
     &                      V_ADV_LBC, V_ADV_LBC_TEND,                  &
     &                      W_ADV_LBC, W_ADV_LBC_TEND,                  &
     &                      MURK_LBC, MURK_LBC_TEND,                    &
     &                      TRACER_LBC, TRACER_LBC_TEND,                &
     &                      1, 0, ErrorStatus, CMESSAGE)

         End If ! L_update_lbcs

        IF ( L_lbc_new ) Then

        !--------------------------------------------------------------
        !           Update primary fields with LAM LBC data
        !--------------------------------------------------------------

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UPDATE_LAM_LBCS',3)
! DEPENDS ON: update_lam_lbcs
        CALL UPDATE_LAM_LBCS(                                                &
     &    r_rho_levels, r_theta_levels,                                 &
     &    ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,               &
     &    OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY,                         &
     &    L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,                 &
     &    L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc, &
     &    L_murk, L_murk_lbc,                                           &
     &    L_LBC_balance, L_int_uvw_lbc,                                 &
     &    RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                        &
     &    LENRIMA(1,1,rima_type_norm),                                  &
     &    LBC_SIZEA(1,1,1,rima_type_norm),                              &
     &    LBC_STARTA(1,1,1,rima_type_norm),                             &
     &    THETA_LBC, Q_LBC, QCL_LBC,                                    &
     &    QCF_LBC, QCF2_LBC, QRAIN_LBC,                                 &
     &    QGRAUP_LBC, CF_BULK_LBC, CF_LIQUID_LBC,                       &
     &    CF_FROZEN_LBC, RHO_LBC,EXNER_LBC,                             &
     &    U_LBC, V_LBC, W_LBC,                                          &
     &    U_ADV_LBC, V_ADV_LBC, W_ADV_LBC,                              &
     &    MURK_LBC,                                                     &
     &    THETA, Q,QCL, QCF,                                            &
     &    QCF2, QRAIN, QGRAUP,                                          &
     &    CF_BULK, CF_LIQUID, CF_FROZEN,                                &
     &    RHO, EXNER_RHO_LEVELS,                                        &
     &    U, V, W, U_ADV, V_ADV, W_ADV, MURK,                           &
     &    DELTA_PHI, DELTA_LAMBDA,                                      &
     &    BASE_PHI, BASE_LAMBDA,                                        &
     &    TWO_OMEGA, DATASTART,                                         &
     &    lat_rot_NP,                                                   &
     &    GLOBAL_ROW_LENGTH, GLOBAL_ROWS                                &
     &     )

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UPDATE_LAM_LBCS',4)

        END IF ! L_lbc_new

      ENDIF     !   model_domain  ==  mt_lam

! Clear workspace
        DEALLOCATE (dtheta_dr_term)

      End If !  ErrorStatus == 0

      Enddo ! end iterations for trajectory calc

      If ( ErrorStatus == 0 ) Then

        Deallocate( u_np1 )
        Deallocate( v_np1 )
        Deallocate( w_np1 )
        Deallocate( theta_np1 )
        Deallocate( rho_np1 )
        Deallocate( q_np1 )
        Deallocate( qcl_np1 )
        Deallocate( qcf_np1 )
        Deallocate( qcf2_np1 )
        Deallocate( qrain_np1 )
        Deallocate( qgraup_np1 )
        Deallocate( mix_v_np1 )
        Deallocate( mix_cl_np1 )
        Deallocate( mix_cf_np1 )
        Deallocate( mix_cf2_np1 )
        Deallocate( mix_rain_np1 )
        Deallocate( mix_graup_np1 )
        Deallocate( mix_v_star)
        Deallocate( mix_cl_star)
        Deallocate( mix_cf_star)
        Deallocate( mix_cf2_star)
        Deallocate( mix_rain_star)
        Deallocate( mix_graup_star)
        If ( NumCycles > 1 ) Then
          Deallocate( R_u_phys1 )
          Deallocate( R_v_phys1 )
          Deallocate( thetastar_phys1 )
          Deallocate( qstar_phys1 )
          Deallocate( qclstar_phys1 )
          Deallocate( qcfstar_phys1 )
          If (L_mcr_qcf2)                                               &
     &      Deallocate( qcf2_star_phys1 )
          If (L_mcr_qrain)                                              &
     &      Deallocate( qrain_star_phys1 )
          If (L_mcr_qgraup)                                             &
     &      Deallocate( qgraup_star_phys1 )
          Deallocate( bulk_cld_frac_phys1 )
          Deallocate( bulk_cld_liq_phys1 )
          Deallocate( bulk_cld_fr_phys1 )
          Deallocate( area_cld_frac_phys1 )
          Deallocate( ti_phys1 )
          Deallocate( zh_phys1 )
          Deallocate( z0msea_phys1 )
          Deallocate( cca_phys1 )
          Deallocate( ccb_phys1 )
          Deallocate( cct_phys1 )
          If ( L_ctile ) Then
            Deallocate( T_LAND_CTILE_PHYS1 )
            Deallocate( T_SICE_CTILE_PHYS1 )
          End If
          Deallocate( t_surf_phys1 )
          Deallocate( t_sf_tile_phys1 )
          Deallocate( snow_tile_phys1 )
          Deallocate( dolr_phys1 )

        End If

!     Deallocate extra microphysics variables
        DEALLOCATE (qcf2_star)
        DEALLOCATE (qrain_star)
        DEALLOCATE (qgraup_star)
        DEALLOCATE (exner_prime)

        If (L_moist_nonhydro_conserve)then
          IF(l_mix_ratio)then
            allocate ( mix_v (1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_levels) )
            allocate ( mix_cl(1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_levels) )
            allocate ( mix_cf(1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_levels) )
            if(L_mcr_qcf2)then
              allocate ( mix_cf2(1-halo_i:row_length+halo_i,            &
     &                           1-halo_j:rows+halo_j, wet_levels) )
            else
              allocate ( mix_cf2(1,1,1) )
            endif
            if(L_mcr_qrain)then
              allocate ( mix_rain(1-halo_i:row_length+halo_i,           &
     &                            1-halo_j:rows+halo_j, wet_levels) )
            else
              allocate ( mix_rain(1,1,1) )
            endif
            if(L_mcr_qgraup)then
              allocate ( mix_graup(1-halo_i:row_length+halo_i,          &
     &                             1-halo_j:rows+halo_j, wet_levels) )
            else
              allocate ( mix_graup(1,1,1) )
            endif
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_v(i,j,k) = mix_v_phys1(i,j,k) + mix_v_phys2(i,j,k)
               mix_cl(i,j,k) = mix_cl_phys1(i,j,k) + mix_cl_phys2(i,j,k)
               mix_cf(i,j,k) = mix_cf_phys1(i,j,k) + mix_cf_phys2(i,j,k)
                End Do
              End Do
            End Do
            if(L_mcr_qcf2)then
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    mix_cf2(i,j,k) = mix_cf2_phys1(i,j,k) +             &
     &                               mix_cf2_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            if(L_mcr_qrain)then
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    mix_rain(i,j,k) = mix_rain_phys1(i,j,k) +           &
     &                                mix_rain_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            if(L_mcr_qgraup)then
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    mix_graup(i,j,k) = mix_graup_phys1(i,j,k) +         &
     &                                 mix_graup_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
! DEPENDS ON: mix_to_q
            call mix_to_q(                                              &
     &                   row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j,                                &
     &                   mix_v, mix_cl, mix_cf,                         &
     &                   mix_cf2, mix_rain, mix_graup,                  &
!    &                   .false.,.false.,.false.,
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   Q, QCL, QCF,                                   &
     &                   QCF2, QRAIN, QGRAUP  )

            if(L_pc2)then
              Do k = 1, wet_levels
                Do j = j_start, j_stop
                  Do i = i_start, i_stop
                    cf_bulk(i,j,k)   = cf_phys1(i,j,k) +cf_phys2(i,j,k)
                    cf_liquid(i,j,k) = cfl_phys1(i,j,k)+cfl_phys2(i,j,k)
                    cf_frozen(i,j,k) = cff_phys1(i,j,k)+cff_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif

            DEALLOCATE (mix_v_phys1)
            DEALLOCATE (mix_cl_phys1)
            DEALLOCATE (mix_cf_phys1)
            DEALLOCATE (mix_v_phys2)
            DEALLOCATE (mix_cl_phys2)
            DEALLOCATE (mix_cf_phys2)
            DEALLOCATE (mix_v)
            DEALLOCATE (mix_cl)
            DEALLOCATE (mix_cf)
            DEALLOCATE (cf_phys1)
            DEALLOCATE (cfl_phys1)
            DEALLOCATE (cff_phys1)
            DEALLOCATE (cf_phys2)
            DEALLOCATE (cfl_phys2)
            DEALLOCATE (cff_phys2)
            DEALLOCATE (mix_cf2_phys1)
            DEALLOCATE (mix_rain_phys1)
            DEALLOCATE (mix_graup_phys1)
            DEALLOCATE (mix_cf2_phys2)
            DEALLOCATE (mix_rain_phys2)
            DEALLOCATE (mix_graup_phys2)
            DEALLOCATE (mix_cf2)
            DEALLOCATE (mix_rain)
            DEALLOCATE (mix_graup)
          else
            Do k = 1, wet_levels
              Do j = j_start, j_stop
                Do i = i_start, i_stop
                  q(i,j,k) = q_phys1(i,j,k) + q_phys2(i,j,k)
                End Do
              End Do
            End Do

            Do k = 1, wet_levels
              Do j = j_start, j_stop
!CDIR NODEP
                Do i = i_start, i_stop
                 qcl(i,j,k) = qcl_phys1(i,j,k) + qcl_phys2(i,j,k)
                 qcf(i,j,k) = qcf_phys1(i,j,k) + qcf_phys2(i,j,k)
                End Do
              End Do
            End Do
            if(L_mcr_qcf2)then
              Do k = 1, wet_levels
                Do j = j_start, j_stop
                  Do i = i_start, i_stop
                    qcf2(i,j,k) = qcf2_phys1(i,j,k) + qcf2_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            if(L_mcr_qrain)then
              Do k = 1, wet_levels
                Do j = j_start, j_stop
                  Do i = i_start, i_stop
                  qrain(i,j,k) = qrain_phys1(i,j,k) + qrain_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            if(L_mcr_qgraup)then
              Do k = 1, wet_levels
                Do j = j_start, j_stop
                  Do i = i_start, i_stop
                qgraup(i,j,k) = qgraup_phys1(i,j,k) + qgraup_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            if(L_pc2)then
              Do k = 1, wet_levels
                Do j = j_start, j_stop
                  Do i = i_start, i_stop
                    cf_bulk(i,j,k)   = cf_phys1(i,j,k) +cf_phys2(i,j,k)
                    cf_liquid(i,j,k) = cfl_phys1(i,j,k)+cfl_phys2(i,j,k)
                    cf_frozen(i,j,k) = cff_phys1(i,j,k)+cff_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            DEALLOCATE (cf_phys1)
            DEALLOCATE (cfl_phys1)
            DEALLOCATE (cff_phys1)
            DEALLOCATE (cf_phys2)
            DEALLOCATE (cfl_phys2)
            DEALLOCATE (cff_phys2)
            DEALLOCATE (q_phys1)
            DEALLOCATE (qcl_phys1)
            DEALLOCATE (qcf_phys1)
            DEALLOCATE (q_phys2)
            DEALLOCATE (qcl_phys2)
            DEALLOCATE (qcf_phys2)
            DEALLOCATE (qcf2_phys1)
            DEALLOCATE (qcf2_phys2)
            DEALLOCATE (qrain_phys1)
            DEALLOCATE (qrain_phys2)
            DEALLOCATE (qgraup_phys1)
            DEALLOCATE (qgraup_phys2)
          endif    !L_mix_ratio

! Use q_pos to remove any negative qcf values
! qlimit being 0.0 try to do locally

! DEPENDS ON: q_pos_ctl
          call Q_Pos_Ctl(                                               &
     &                   QCF, row_length, rows, wet_levels,             &
     &                   global_row_length, global_rows,                &
     &                   mype, nproc, halo_i, halo_j,                   &
     &                   gc_all_proc_group,                             &
     &                   model_domain,                                  &
     &                   halo_type_extended, l_q_pos_local, 0.0         &
     &                   )

        Else If ( NumCycles >1 .and. L_pc2 ) Then
! If not L_moist_nonhydro_conserve but cycling and L_pc2
! is used deallocate cf_phys1, cfl_phys1, cff_phys1
          DEALLOCATE (cf_phys1)
          DEALLOCATE (cfl_phys1)
          DEALLOCATE (cff_phys1)
        End If ! L_moist_nonhydro_conserve



      End If        !  ErrorStatus  ==  0

! ----------------------------------------------------------------------
! section 6.2  Check for q below qlimit and reset
! ----------------------------------------------------------------------

      If(L_qpos)then

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('Q_Pos_Ctl',3)
! DEPENDS ON: q_pos_ctl
        call Q_Pos_Ctl(                                                 &
     &                   Q, row_length, rows, wet_levels,               &
     &                   global_row_length, global_rows,                &
     &                   mype, nproc, halo_i, halo_j,                   &
     &                   gc_all_proc_group,                             &
     &                   model_domain,                                  &
     &                   halo_type_extended, l_q_pos_local, qlimit      &
     &                   )
!
        If (L_pc2) then
! For the PC2 cloud scheme we also need to limit the liquid

! DEPENDS ON: q_pos_ctl
          call Q_Pos_Ctl(                                               &
     &                   QCL, row_length, rows, wet_levels,             &
     &                   global_row_length, global_rows,                &
     &                   mype, nproc, halo_i, halo_j,                   &
     &                   gc_all_proc_group,                             &
     &                   model_domain,                                  &
     &                   halo_type_extended, l_q_pos_local, 0.0         &
     &                   )
!
       End If       ! L_pc2
!
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('Q_Pos_Ctl',4)

      End If       ! L_qpos

      if(sf(0,30)) then
        if(sf(182,30).or.sf(183,30).or.sf(184,30).or.sf(189,30).or.     &
     &     sf(190,30).or.sf(191,30) ) then
! increment calculations for moisture variables
          Do k = 1, wet_levels
            Do j = j_start, j_stop
              Do i = i_start, i_stop
               inc_q(i,j,k)   = q(i,j,k)   - inc_q(i,j,k)
               inc_qcl(i,j,k) = qcl(i,j,k) - inc_qcl(i,j,k)
               inc_qcf(i,j,k) = qcf(i,j,k) - inc_qcf(i,j,k)
               inc_cf(i,j,k)  = cf_bulk(i,j,k)   - inc_cf(i,j,k)
               inc_cfl(i,j,k) = cf_liquid(i,j,k) - inc_cfl(i,j,k)
               inc_cff(i,j,k) = cf_frozen(i,j,k) - inc_cff(i,j,k)
              End Do
            End Do
          End Do
        endif
      endif      !    sf(0,30)

! ----------------------------------------------------------------------
! Section 7.0 Mean all polar variables on a level to remove deviation
!             of values due to rounding error.
!             Only called every polar_reset_timesteps
! ----------------------------------------------------------------------
      If (ErrorStatus  ==  0 .and. model_domain  ==  mt_global          &
     &    .and. L_polar_reset) Then

! need if test on correct timestep
        If (Timestep_Number  ==  1 .or. (                               &
     &      (Timestep_number-1)/polar_reset_timesteps*                  &
     &      polar_reset_timesteps  ==  Timestep_number-1 ) )            &
     &    Then

! DEPENDS ON: polar_reset_mean
          Call Polar_Reset_Mean(                                        &
     &                      EXNER_RHO_LEVELS,RHO,                       &
     &                      THETA,W,                                    &
     &                      Q, QCL, QCF,                                &
     &                      CF_BULK, CF_LIQUID,                         &
     &                      CF_FROZEN,                                  &
     &                      row_length, rows, model_levels,             &
     &                      wet_levels, global_row_length,              &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      nproc, nproc_y, gc_proc_row_group,          &
     &                      at_extremity)
        End If
      End If

! ------------------------------------------------------------------
! Section 17  Aerosol Modelling - includes Sulphur Cycle and Soot
!
! ------------------------------------------------------------------
!
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AEROSOL MODELLING',5)
!
      IF (L_SULPC_SO2.OR.L_SOOT.OR.L_BIOMASS.OR.L_OCFF) THEN
!
! Allocate diagnostic space for STASH
        ALLOCATE (STASHwork17(STASH_maxlen(17,A_im)))
!
! Don't call Swap_bounds for fields used in Aero_Ctl
!
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('AEROSOL MODELLING',3)
!
        ALLOCATE(O3_MMR(row_length,rows,model_levels))
        ALLOCATE(H2O2_MMR(row_length,rows,model_levels))
        ALLOCATE(OH_conc(row_length,rows,model_levels))
        ALLOCATE(HO2_conc(row_length,rows,model_levels)) 
!

! DEPENDS ON: get_sulpc_oxidants
        CALL GET_SULPC_OXIDANTS(                                        &
     & L_SULPC_ONLINE_OXIDANTS, L_UKCA, L_UKCA_TROP, L_UKCA_TROPISOP,   &
     & L_UKCA_STRATTROP, first_atmstep_call,                            &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Start arg_atm_fields.h

! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!
! Current Code Owner: A. Treshansky
!


      ! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP, &
      ! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV, &
      ! 1.2: Data variables stored in secondary space.
       P, & 
      ! Pressure on theta levels
       P_THETA_LEVELS, &
      ! Exner pressure on theta levels
       EXNER_THETA_LEVELS, &
      ! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN, &
      ! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF, &
      ! 1.5: Radiation Increments
       SW_INCS, LW_INCS, &
! PAR radiation increment
       DIRPAR, &
      ! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP, &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3, &
!  tropopause-based ozone
       TPPSOZONE, &
      ! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK, &       
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6, &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3, &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD, &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, &
       CO2, CH4_STOCH, O3_STOCH, &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5, &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10, &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15, &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20, &
! CABLE
!       TSOIL_TILE, SMCL_TILE, STHU_TILE, STHF_TILE, SNOW_DEPTH3L,       &
       TSOIL_TILE, SMCL_TILE,  STHF_TILE, SNOW_DEPTH3L,                 &
       SNOW_MASS3L, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SNOW_AGE,       & 
       SNOW_FLG3L,                                                      &
! Lestevens Sept 2012 - adding progs for CASACNP
       CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
       NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
       cable_lai,                                                       &

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC, &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC, &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC, &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC, TRACER_LBC, &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND, &
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND, &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND, &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND, &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND, &
       MURK_LBC_TEND, TRACER_LBC_TEND, &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB, &
      ! 2.2: Data variables stored in secondary space.
       PSTAR, &
      ! 2.3: Cloud fields
       CCB, CCT, CCLWP, &
      ! 2.4: Boundary layer fields
       ZH, & 
      ! Standard deviation of turbulent fluctuations of layer 1
       T1_SD, &
      ! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD, &
      ! Number of model levels in the  turbulently mixed layer
       NTML, &
      ! Top level for turb mixing in any decoupled Sc layer
       NTDSC, &
      ! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS, & 
      ! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT, &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN, &
      ! 2.5: Vegetation Ancillary fields
       CANOPY_WATER, SURF_CAP, SURF_RESIST, ROOT_DEPTH, INFILT, &
       VEG_FRAC, LAI, CANHT, Z0, SFA, MDSA, GS, &
      ! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2, &
       OROG_GRAD_X, OROG_GRAD_Y, &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY, &
      ! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS, &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, &
      ! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT, &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND, &
       DUST_MREL1, DUST_MREL2, DUST_MREL3, &
       DUST_MREL4, DUST_MREL5, DUST_MREL6, &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM, &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX, &
      ! Tracer Fluxes - kdcorbin, 05/10
       TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4, &
       TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8, &
       TRACER_FLUX9, TRACER_FLUX10, TRACER_FLUX11, TRACER_FLUX12, &
       TRACER_FLUX13, TRACER_FLUX14, TRACER_FLUX15, TRACER_FLUX16, &
       TRACER_FLUX17, TRACER_FLUX18, TRACER_FLUX19, TRACER_FLUX20, &

      ! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5, &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10, &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15, &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20, &
      !   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX, &
      !   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5, &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9, &
       LAI_PFT, CANHT_PFT, DISTURB_VEG, &
       SOIL_ALB, SOIL_CARB, &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4, &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC, &
       RSP_W_PFT_ACC, RSP_S_ACC, &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4, &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE, &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, &
       DOLR_FIELD, &
       LW_DOWN, SW_TILE_RTS, &
!! MODIFIED BY AT
!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!       TSLAB, TCLIM, HCLIM, CHEAT, OIFLX, UICE, VICE, &
!       SIG11NE, SIG11SE, SIG11SW, SIG11NW, &
!       SIG12NE, SIG12SE, SIG12SW, SIG12NW, &
!       SIG22NE, SIG22SE, SIG22SW, SIG22NW, &
!! END MODIFIED BY AT
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS, &
!   2.15: Fields carried forward from previous version.
!         May not be required
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho, &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP, &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND, &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET, &
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE, &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM, &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT, &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN, &
       C_SOLAR,C_BLUE,C_DOWN,C_LONGWAVE,C_TAUX,C_TAUY,C_WINDMIX, &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_BOTMELTN,C_TOPMELTN,C_LSRAIN, &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_PRESS, C_U10, C_V10, &
! UKCA oxidant fields
        OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA, & 
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR, &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK, &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4, &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC, &
       ARCLDLTA_DL, &
! Convective Cloud Fields
       LCBASE, CCW_RAD, &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM, &

! End arg_atm_fields.h
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
     & O3_MMR, H2O2_MMR, OH_conc, HO2_conc)
!
! DEPENDS ON: aero_ctl
        CALL AERO_CTL(                                                  &
! Parallel variables
     &  halo_i, halo_j, offx, offy, global_row_length, global_rows      &
     &, gc_proc_row_group, gc_proc_col_group                            &
     &, at_extremity, nproc, nproc_x, nproc_y                           &
     &, neighbour, g_rows, g_row_length, g_datastart, mype              &
! model dimensions
     &, row_length, rows, n_rows, land_field                            &
     &, model_levels, wet_levels, bl_levels, n_cca_lev                  &
     &, theta_field_size                                                &
     &, salt_dim1, salt_dim2, salt_dim3                                 &
     &, aero_dim1, aero_dim2, aero_dim3                                 &
! Model switches
     &, model_domain, LCAL360, L_SEC_VAR, L_EqT, Ltimer                 &
! Model parameters
     &, Ntot_land, Ntot_sea                                             &
! Physical constants
     &, Lc, Lf, Cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse, earth_radius, Pi                                   &
! Co-ordinate information
     &, r_rho_levels, r_theta_levels                                    &
     &, eta_theta_levels, eta_rho_levels                                &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &
! Time stepping information
     &, timestep                                                        &
     &, I_year, I_day_number, I_hour, I_minute                          &
     &, I_second, timestep_number                                       &
     &, PREVIOUS_TIME                                                   &
     &, CALL_CHEM_FREQ                                                  &
! Trig arrays
     &, sin_theta_longitude, cos_theta_longitude                        &
     &, FV_cos_theta_latitude                                           &
! Grid-dependent arrays
     &,     f3_at_u, true_longitude, true_latitude                      &
!
! Data fields IN
     &, U, V, TSTAR, TSTAR_SEA                                          &
     &, THETA, Q, QCL, QCF                                              &
     &, RHO, LAND, FRAC_LAND, P_THETA_LEVELS                            &
     &, EXNER_RHO_LEVELS, EXNER_THETA_LEVELS                            &
     &, ICE_FRACTION, SNODEP                                            &
     &, CF_BULK                                                         &
     &, OH_conc, H2O2_MMR, HO2_conc, O3_MMR                             &
     &, SO2_EM, SO2_HILEM, SO2_NATEM                                    &
     &, DMS_EM, DMS_CONC, NH3_EM                                        &
     &, DMS_OFLUX                                                       &
     &, SOOT_EM, SOOT_HILEM, BMASS_EM, BMASS_HILEM, OCFF_EM, OCFF_HILEM &
     &, SO2_HIGH_LEVEL, SOOT_HIGH_LEVEL, BMASS_HIGH_LEVEL_1             &
     &, BMASS_HIGH_LEVEL_2, OCFF_HIGH_LEVEL, land_index                 &
! Logicals IN
     &, L_SULPC_SO2, L_SULPC_DMS, L_SULPC_OZONE                         &
     &, L_SULPC_SO2_O3_NONBUFFERED, L_SULPC_NH3                         &
     &, L_use_sulphate_sulpc, L_use_seasalt_sulpc, L_SOOT               &
     &, l_use_soot_sulpc, l_biomass, l_use_bmass_sulpc                  &
     &, l_ocff, l_use_ocff_sulpc                                        &
     &, L_SO2_SURFEM, L_SO2_HILEM, L_SO2_NATEM, L_DMS_EM                &
     &, L_DMS_em_inter, L_DMS_Liss_Merlivat                             &
     &, L_DMS_Wanninkhof, L_DMS_Nightingale                             &
     &, L_DMS_Ointer                                                    &
     &, L_NH3_EM, L_CTILE                                               &
     &, L_SOOT_SUREM, L_SOOT_HILEM, L_BMASS_SUREM, L_BMASS_HILEM        &
     &, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_BIOGENIC                      &
     &, L_USE_SEASALT_DIRECT, L_USE_SEASALT_INDIRECT                    &
     &, L_USE_SEASALT_AUTOCONV, L_DUST                                  &
!
! Data fields IN/OUT
     &, SO2, DMS                                                        &
     &, SO4_AITKEN, SO4_ACCU, SO4_DISS                                  &
     &, H2O2, NH3                                                       &
     &, SOOT_NEW, SOOT_AGD, SOOT_CLD                                    &
     &, BMASS_NEW, BMASS_AGD, BMASS_CLD                                 &
     &, OCFF_NEW, OCFF_AGD, OCFF_CLD                                    &
     &, biogenic                                                        &
!
! Data fields IN
     &, DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5           &
!
! Data fields OUT
! Diagnostic info
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &  STASHwork17                                                     &
! Error info
     &, ErrorStatus                                                     &
     & )
!
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('AEROSOL MODELLING',4)
!
! Don't call Swap_bounds for updated fields
!
! Diagnostics STASHed for Aerosol section 17
!
! DEPENDS ON: timer
        If (Ltimer) CALL TIMER('STASH',3)
!
! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,17,STASHwork17,                            &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    ErrorStatus,Cmessage)
!
! DEPENDS ON: timer
        If (Ltimer) CALL TIMER('STASH',4)
!
        DEALLOCATE (STASHwork17)
!
        IF(ALLOCATED(O3_MMR)) DEALLOCATE(O3_MMR)
        IF(ALLOCATED(H2O2_MMR)) DEALLOCATE(H2O2_MMR)
        IF(ALLOCATED(OH_conc)) DEALLOCATE(OH_conc)
        IF(ALLOCATED(HO2_conc)) DEALLOCATE(HO2_conc)

      END IF         ! END L_SULPC_SO2.OR.L_SOOT TEST
!
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AEROSOL MODELLING',6)
!
! ----------------------------------------------------------------------
! Section 8.0 Update pressure to value at new-time level.
!             Calculate surface pressure, and exner and pressure at
!             theta levels.
!             If old lbcs then need to update exner on boundaries first.
! ----------------------------------------------------------------------
      If (ErrorStatus  ==  0 .and. model_domain  ==  mt_lam) Then

       If ( .NOT. L_lbc_new ) Then

! Calculate EXNER_LBCs at next time level

! Change halo size for exner_lbc
        lbc_size=LENRIMA(fld_type_p,halo_type_extended,rima_type_norm)

        If (RIM_STEPSA  ==  0) Then

          ! No LBC updating
          Do k=1,MODEL_LEVELS+1
            Do i=1,lbc_size
              EXNER_LBC_REAL_TEND(i,k)=EXNER_LBC(i,k)  
            End Do !i
          End Do !k

        Else If (MOD(Timestep_Number,RIM_STEPSA)  ==  0) Then

          ! End of current LBC period
          DO k=1,MODEL_LEVELS+1
            DO i=1,lbc_size
              EXNER_LBC_REAL_TEND(i,k)= EXNER_LBC_TEND(i,k)
            ENDDO !i
          ENDDO !k

        ELSE ! Just a normal timestep during a LBC period

          increment_factor=1.0/                                         &
     &      (RIM_STEPSA-MOD(Timestep_Number-1,RIM_STEPSA))

          DO k=1,MODEL_LEVELS+1
            DO i=1,lbc_size
              EXNER_LBC_REAL_TEND(i,k) = EXNER_LBC(i,k) +               &
     &                                              increment_factor *  &
     &                            (EXNER_LBC_TEND(i,k) - EXNER_LBC(i,k))
            ENDDO !i
          ENDDO !k

        ENDIF ! End of current LBC period?

!     n_rims_to_do set at begining of subroutine. Only wts=1 rims
!     updated here since derived p fields are only being used for
!     diagnostics. They get properly weighted at start of next timestep
        L_do_halos=.FALSE.
        L_do_boundaries=.TRUE.

! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
     &    ROW_LENGTH,ROWS,Offx,Offy,                                    &
     &    MODEL_LEVELS+1,fld_type_p,EXNER_RHO_LEVELS,                   &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    EXNER_LBC_REAL_TEND,                                          &
     &    RIMWIDTHA(rima_type_norm),n_rims_to_do,                       &
     &    RIMWEIGHTSA,AT_EXTREMITY,                                     &
     &    L_do_boundaries,L_do_halos)

       End If ! .NOT. L_lbc_new

      Endif      ! ErrorStatus  ==  0 .and. model_domain  ==  mt_lam

      If (ErrorStatus  ==  0) Then

! Check for negative pressure if requested to
        If (Instability_diagnostics  >   0) Then
          Do k = 1, model_levels+1
            Do j = 1, rows
              Do i = 1, row_length
                ij = i+offx + (j+offy-1) * (row_length+2*offx)
                If (exner_rho_levels(((k-1)*theta_off_size)+ij) <  0.) Then
                  ErrorStatus = 123
                End If
              End Do
            End Do
          End Do
! ErrorStatus 123 message now generated outside of the loops.
          if ( ErrorStatus == 123 ) THEN
! DEPENDS ON: ereport
            Call Ereport("ATM_STEP", ErrorStatus,                       &
     &           "Negative pressure value about to be created" )
          endif
        End If       !    Instability_diagnostics  >   0

      End If      ! ErrorStatus  ==  0

      If ( ErrorStatus  ==  0) Then
! DEPENDS ON: consistent_pressure
        call     Consistent_Pressure (                                  &
     &           exner_rho_levels,                                      &
     &           offx,offy,halo_i,halo_J,                               &
     &           row_length,rows,model_levels,                          &
     &           kappa, g, r_theta_levels, r_rho_levels, rho, p_zero,   &
     &           p, pstar, p_theta_levels,exner_theta_levels)
!
        IF (L_physics) then
!
! Are we using the PC2 cloud scheme?
!
        If (L_pc2) then
!
! ----------------------------------------------------------------------
! PC2: Calculate condensation due to changes in temperature resulting
!    from adiabatic changes in pressure (mainly from vertical advection)
!    Also call checking routine from within this subroutine.
! ----------------------------------------------------------------------
!
          allocate(t_inc_pres  (row_length,rows,wet_levels))
          allocate(q_inc_pres  (row_length,rows,wet_levels))
          allocate(qcl_inc_pres(row_length,rows,wet_levels))
          allocate(qcf_inc_pres(row_length,rows,wet_levels))
          allocate(cf_inc_pres (row_length,rows,wet_levels))
          allocate(cfl_inc_pres(row_length,rows,wet_levels))
          allocate(cff_inc_pres(row_length,rows,wet_levels))
          allocate(t_dini  (row_length,rows,wet_levels))
          allocate(q_dini  (row_length,rows,wet_levels))
          allocate(qcl_dini(row_length,rows,wet_levels))
          allocate(qcf_dini(row_length,rows,wet_levels))
          allocate(cf_dini (row_length,rows,wet_levels))
          allocate(cfl_dini(row_length,rows,wet_levels))
          allocate(cff_dini(row_length,rows,wet_levels))

! DEPENDS ON: pc2_pressure_forcing
          Call pc2_pressure_forcing(                                    &
     &              halo_i, halo_j, offx, offy,                         &
     &              P,PSTAR,                                            &
     &              P_THETA_LEVELS, wet_levels,                         &
     &              row_length, rows, rhc_row_length, rhc_rows,         &
     &              timestep, rhcpt, THETA, CF_BULK,                    &
     &              CF_LIQUID, CF_FROZEN,                               &
     &              Q, QCL,QCF,                                         &
     &              exner_star, EXNER_THETA_LEVELS,                     &
     &              CCB, CUMULUS, rhts, tlts, qtts, ptts, CF_AREA,      &
     &              t_inc_pres, q_inc_pres, qcl_inc_pres, qcf_inc_pres, &
     &              cf_inc_pres, cfl_inc_pres,                          &
     &   cff_inc_pres, t_dini, q_dini, qcl_dini, qcf_dini,              &
     &   cf_dini, cfl_dini, cff_dini, l_mr_pc2, L_ACF_Cusack )

          deallocate(rhts)
          deallocate(tlts)
          deallocate(qtts)
          deallocate(ptts)

        End if  ! L_pc2

          If ( (.not. L_pc2) .or. l_pc2_reset ) Then
! ----------------------------------------------------------------------
! Call cloud scheme to make cloud consistent with moisture fields
! ----------------------------------------------------------------------

! DEPENDS ON: qt_bal_cld
            call qt_bal_cld(                                            &
     &        PSTAR,P_THETA_LEVELS,P,                                   &
     &        THETA,EXNER_THETA_LEVELS,                                 &
     &        Q,QCL,QCF,QCF2,                                           &
     &        rhcpt, rhc_row_length, rhc_rows, bl_levels,               &
     &        cloud_fraction_method,overlap_ice_liquid,                 &
     &        ice_fraction_method,ctt_weight,t_weight,                  &
     &        qsat_fixed,sub_cld,                                       &
     &        row_length,rows,model_levels,wet_levels,                  &
     &        offx,offy,halo_i,halo_j,                                  &
     &        delta_lambda, delta_phi,                                  &
     &        r_theta_levels, FV_cos_theta_latitude,                    &
     &        lc, cp, L_cld_area, L_ACF_Cusack, L_ACF_Brooks,           &
     &        L_eacf, L_mcr_qcf2,                                       &
     &        l_mr_qtbalcld,                                            &
     &        ntml, cumulus, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,    &
     &        mype)

          Endif   ! L_pc2 and L_pc2_reset
          If (L_run_with_physics2) DEALLOCATE ( RHCPT )
        Endif   ! L_physics

        if (sf(181,30) ) then
! Store increment to T
          Do k = 1, model_levels
            Do j = j_start, j_stop
              Do i = i_start, i_stop
                inc_t(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k) &
     &            - inc_t(i,j,k)
              End Do
            End Do
          End Do
        endif     ! stashflag

      End If      !  ErrorStatus  ==  0

! DEPENDS ON: timer
      If (Ltimer) CALL TIMER('AS5-8 Updates',6)

! ----------------------------------------------------------------------
! Section 9.0a Calculate Mass and energy of atmosphere if required
!              using end of timestep values
!             Only done every energy correction period.
! ----------------------------------------------------------------------
      IF (L_EMCORR) THEN

! Set energy correction for use in section 30
        energy_corr_now=A_REALHD(rh_energy_corr)

        IF (LENERGY) THEN
! DEPENDS ON: timer
          IF(LTIMER) CALL TIMER('AS9 Energy mass   ',5)
! zero fields to be calculated

          tot_energy_final = 0.0
          tot_dry_mass_final = 0.0
          tot_moist_final = 0.0

! DEPENDS ON: eng_mass_diag
          Call eng_mass_diag (                                          &
! Parallel variables
     &                      halo_i, halo_j, offx, offy                  &
     &,                     global_row_length, gc_proc_row_group        &
     &,                     gc_proc_col_group                           &
     &,                     at_extremity, nproc, nproc_x, nproc_y       &
     &,                     neighbour                                   &
     &,                     mype                                        &
! model info
     &,                     row_length, rows, n_rows                    &
     &,                     model_domain                                &
     &,                     model_levels,wet_levels                     &
     &,                     r_theta_levels,r_rho_levels                 &
     &,                     delta_lambda,delta_phi                      &
! trig
     &,                     FV_cos_theta_latitude,cos_v_latitude        &
     &,                     cos_theta_longitude,sin_theta_longitude     &
! data arrays in only
     &,                     THETA , U, V                                &
     &,                     W, RHO , Q                                  &
     &,                     QCL, QCF                                    &
     &,                     EXNER_THETA_LEVELS                          &
!     &,                     d1(jexner_rho_levels(1))
!     &,                     d1(jp(1)), d1(jp_theta_levels(1))
!     &,                     d1(jpstar)
! sum of moist fluxes
     &,                     NET_MFLUX                                   &
     &,                     A_REALHD(rh_tot_mass_init)                  &
     &,                     A_REALHD(rh_tot_m_init)                     &
! logical to indicate mass and moist correction required
     &,                     Lmass_corr,Lqt_corr,Lemq_print              &
! energy correction timestep info
     &,                     a_energysteps,timestep                      &
! IN/OUT  results from calculations
     &,                     tot_energy_final,tot_dry_mass_final         &
     &,                     tot_moist_final)

! DEPENDS ON: timer
          If (Ltimer) Call timer ('cal_eng_mass_corr',3)

! DEPENDS ON: cal_eng_mass_corr
          Call cal_eng_mass_corr (                                      &
     &                      row_length, rows, mype                      &
     &,                     a_energysteps,timestep                      &
     &,                     delta_lambda,delta_phi                      &
     &,                     halo_i, halo_j                              &
     &,                     NET_FLUX                                    &
     &,                     A_REALHD(rh_tot_mass_init)                  &
     &,                     A_REALHD(rh_tot_energy_init)                &
     &,                     A_REALHD(rh_energy_corr)                    &
     &,                     tot_energy_final )

! DEPENDS ON: timer
          If (Ltimer) Call timer ('cal_eng_mass_corr',4)

! Swap initial energy and final energy.

          A_REALHD(rh_tot_energy_init) = tot_energy_final

! Swap initial moisture and final moisture.

          A_REALHD(rh_tot_m_init) = tot_moist_final

! DEPENDS ON: timer
          IF(LTIMER) CALL TIMER('AS9 Energy mass   ',6)

        ENDIF   ! LENERGY

      ELSE
! Set energy correction for use in section 30
        energy_corr_now=0.0
      ENDIF     ! L_EMCORR

! ----------------------------------------------------------------------
! Section ?: If required, call temporal filtering control routine.
! ----------------------------------------------------------------------
      Sec = STEPIM(A_IM) * INT(SECS_PER_STEPIM(A_IM))

      L_CallTFilt = ( L_IAU .AND. (Sec >= IAU_StartMin * 60)            &
     &                      .AND. (Sec <= IAU_EndMin   * 60) )          &
     &              .OR.                                                &
     &              ( L_TDF .AND. (Sec >= TDF_StartMin * 60)            &
     &                      .AND. (Sec <= TDF_EndMin   * 60) )

      IF (L_CallTFilt) THEN

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('TFILTCTL',3)

! DEPENDS ON: tfilt_cntl
        CALL TFilt_cntl (                                               &
      U, V, W, U_adv, V_adv, W_adv,                                     &
      Theta, Exner_rho_levels, Rho,                                     &
      Q, Qcl, Qcf, Murk, ozone_tracer,                                  &
      Deep_soil_temp,                                                   &
      P, P_theta_levels, Exner_theta_levels,                            &
      snodep,                                                           &
      cf_area, cf_bulk, cf_liquid, cf_frozen,                           &
      Pstar, Tstar, Tstar_tile,                                         &
! ARGDUMA Dump headers
     &  A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
     &  A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
     &  A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGCONA start
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field sin_u_latitude. J F Thomson.
!  5.0   02/06/99  Insert C-P C-grid constants. M L Gallani.
!  5.3   01/10/01  Add fields for chequerboard radiation. S Cusack
! 6.1  04/08/04  Add separate arrays diff_coeff_u, diff_coeff_v
!                                                     Terry Davies
! 6.2  05/01/06   Add true_latitude. Yongming Tang
! 6.2  14/12/06  Add separate arrays VarRes Array co-ordinates
!                                                       Terry Davies

! argcona.h contained constants for the atmosphere.
! As of vn6.6 these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, AD_MASK_TROP_MOD, ROT_COEFF_MOD

! ARGLNDM Constants for physics routines
     &  land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &                    l_mr_tfiltctl,                                &
                                         ! In. Use mixing ratios
     &                    IAU_lookup,                                   &
                                      ! inout
     &                    D1_IAU_k4,                                    &
                                      ! inout
     &                    D1_TFilt )  ! inout

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('TFILTCTL',4)

      END IF      !  L_CallTFilt
!
! Are we using the PC2 cloud scheme to determine area cloud fraction?
!
        If (L_pc2 .and. .not. L_pc2_reset) then
          If (.not. L_cld_area) then
!
! ----------------------------------------------------------------------
! PC2: Set area cloud fraction to the bulk cloud fraction. Use the
!    D1 arrays directly
! ----------------------------------------------------------------------
!
          Do k = 1, wet_levels
            Do j = j_start, j_stop
              Do i = i_start, i_stop
                CF_AREA(i,j,k) = CF_BULK(i,j,k)
              End do
            End do
          End do

          Else If (L_cld_area) then

            If (L_ACF_Brooks) then
              Allocate ( cf_bulk_nohalo(row_length,rows,wet_levels) )
              Allocate ( cf_liquid_nohalo(row_length,rows,wet_levels) )
              Allocate ( cf_frozen_nohalo(row_length,rows,wet_levels) )

! Place bulk, liquid and frozen cloud fractions in halo-free arrays
! Use indexing over the full row and row_length (including any LAM
! boundary rim) since the call to ls_acf_brooks uses this indexing.
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    cf_bulk_nohalo(i,j,k)   = cf_bulk(i,j,k)
                    cf_liquid_nohalo(i,j,k) = cf_liquid(i,j,k)
                    cf_frozen_nohalo(i,j,k) = cf_frozen(i,j,k)
                  End Do
                End Do
              End Do

! DEPENDS ON: ls_acf_brooks
              Call LS_ACF_Brooks (                                      &
     &             halo_i, halo_j, Offx, Offy                           &
     &,            row_length, rows, model_levels, wet_levels           &
     &,            r_theta_levels, delta_lambda, delta_phi              &
     &,            FV_cos_theta_latitude                                &
     &,            cf_bulk_nohalo, cf_liquid_nohalo                     &
     &,            cf_frozen_nohalo, cumulus                            &
     &,            CF_AREA )

              Deallocate ( cf_bulk_nohalo )
              Deallocate ( cf_liquid_nohalo )
              Deallocate ( cf_frozen_nohalo )

            End If ! L_ACF_Brooks

          End If ! L_cld_area
!
        End if  ! L_pc2 .and. .not. L_pc2_reset

! ----------------------------------------------------------------------
! Section 9.0 Diagnostics at end of timestep
! ----------------------------------------------------------------------

! DEPENDS ON: timer
      IF(LTIMER) CALL TIMER('AS9 End TStep Diags',5)

! section 15: 'dynamics' based quantities
      IF(      SF(0,15)                                                 &
                                ! Diagnostics required for this section
     &   .AND. ErrorStatus == 0) THEN

! DEPENDS ON: st_diag1
        CALL St_diag1(STASH_maxlen(15,A_im),                            &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
     &  A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
     &  A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
     &  A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Start arg_atm_fields.h

! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!
! Current Code Owner: A. Treshansky
!


      ! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP, &
      ! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV, &
      ! 1.2: Data variables stored in secondary space.
       P, & 
      ! Pressure on theta levels
       P_THETA_LEVELS, &
      ! Exner pressure on theta levels
       EXNER_THETA_LEVELS, &
      ! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN, &
      ! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF, &
      ! 1.5: Radiation Increments
       SW_INCS, LW_INCS, &
! PAR radiation increment
       DIRPAR, &
      ! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP, &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3, &
!  tropopause-based ozone
       TPPSOZONE, &
      ! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK, &       
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6, &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3, &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD, &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, &
       CO2, CH4_STOCH, O3_STOCH, &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5, &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10, &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15, &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20, &
! CABLE
!       TSOIL_TILE, SMCL_TILE, STHU_TILE, STHF_TILE, SNOW_DEPTH3L,       &
       TSOIL_TILE, SMCL_TILE,  STHF_TILE, SNOW_DEPTH3L,                 &
       SNOW_MASS3L, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SNOW_AGE,       & 
       SNOW_FLG3L,                                                      &
! Lestevens Sept 2012 - adding progs for CASACNP
       CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
       NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
       cable_lai,                                                       &

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC, &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC, &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC, &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC, TRACER_LBC, &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND, &
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND, &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND, &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND, &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND, &
       MURK_LBC_TEND, TRACER_LBC_TEND, &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB, &
      ! 2.2: Data variables stored in secondary space.
       PSTAR, &
      ! 2.3: Cloud fields
       CCB, CCT, CCLWP, &
      ! 2.4: Boundary layer fields
       ZH, & 
      ! Standard deviation of turbulent fluctuations of layer 1
       T1_SD, &
      ! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD, &
      ! Number of model levels in the  turbulently mixed layer
       NTML, &
      ! Top level for turb mixing in any decoupled Sc layer
       NTDSC, &
      ! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS, & 
      ! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT, &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN, &
      ! 2.5: Vegetation Ancillary fields
       CANOPY_WATER, SURF_CAP, SURF_RESIST, ROOT_DEPTH, INFILT, &
       VEG_FRAC, LAI, CANHT, Z0, SFA, MDSA, GS, &
      ! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2, &
       OROG_GRAD_X, OROG_GRAD_Y, &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY, &
      ! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS, &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, &
      ! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT, &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND, &
       DUST_MREL1, DUST_MREL2, DUST_MREL3, &
       DUST_MREL4, DUST_MREL5, DUST_MREL6, &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM, &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX, &
      ! Tracer Fluxes - kdcorbin, 05/10
       TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4, &
       TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8, &
       TRACER_FLUX9, TRACER_FLUX10, TRACER_FLUX11, TRACER_FLUX12, &
       TRACER_FLUX13, TRACER_FLUX14, TRACER_FLUX15, TRACER_FLUX16, &
       TRACER_FLUX17, TRACER_FLUX18, TRACER_FLUX19, TRACER_FLUX20, &

      ! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5, &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10, &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15, &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20, &
      !   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX, &
      !   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5, &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9, &
       LAI_PFT, CANHT_PFT, DISTURB_VEG, &
       SOIL_ALB, SOIL_CARB, &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4, &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC, &
       RSP_W_PFT_ACC, RSP_S_ACC, &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4, &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE, &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, &
       DOLR_FIELD, &
       LW_DOWN, SW_TILE_RTS, &
!! MODIFIED BY AT
!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!       TSLAB, TCLIM, HCLIM, CHEAT, OIFLX, UICE, VICE, &
!       SIG11NE, SIG11SE, SIG11SW, SIG11NW, &
!       SIG12NE, SIG12SE, SIG12SW, SIG12NW, &
!       SIG22NE, SIG22SE, SIG22SW, SIG22NW, &
!! END MODIFIED BY AT
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS, &
!   2.15: Fields carried forward from previous version.
!         May not be required
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho, &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP, &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND, &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET, &
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE, &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM, &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT, &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN, &
       C_SOLAR,C_BLUE,C_DOWN,C_LONGWAVE,C_TAUX,C_TAUY,C_WINDMIX, &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_BOTMELTN,C_TOPMELTN,C_LSRAIN, &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_PRESS, C_U10, C_V10, &
! UKCA oxidant fields
        OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA, & 
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR, &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK, &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4, &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC, &
       ARCLDLTA_DL, &
! Convective Cloud Fields
       LCBASE, CCW_RAD, &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM, &

! End arg_atm_fields.h
! ARGCONA start
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field sin_u_latitude. J F Thomson.
!  5.0   02/06/99  Insert C-P C-grid constants. M L Gallani.
!  5.3   01/10/01  Add fields for chequerboard radiation. S Cusack
! 6.1  04/08/04  Add separate arrays diff_coeff_u, diff_coeff_v
!                                                     Terry Davies
! 6.2  05/01/06   Add true_latitude. Yongming Tang
! 6.2  14/12/06  Add separate arrays VarRes Array co-ordinates
!                                                       Terry Davies

! argcona.h contained constants for the atmosphere.
! As of vn6.6 these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, AD_MASK_TROP_MOD, ROT_COEFF_MOD

! ARGLNDM Constants for physics routines
     &  land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     & ErrorStatus,CMessage)

      ENDIF ! Diagnostics required for this section

! section 16: 'physics' based quantities
      IF(      SF(0,16)                                                 &
                                ! Diagnostics required for this section
     &   .AND. ErrorStatus == 0) THEN

! DEPENDS ON: st_diag2
        CALL St_diag2(STASH_maxlen(16,A_im),                            &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
     &  A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
     &  A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
     &  A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Start arg_atm_fields.h

! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!
! Current Code Owner: A. Treshansky
!


      ! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP, &
      ! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV, &
      ! 1.2: Data variables stored in secondary space.
       P, & 
      ! Pressure on theta levels
       P_THETA_LEVELS, &
      ! Exner pressure on theta levels
       EXNER_THETA_LEVELS, &
      ! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN, &
      ! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF, &
      ! 1.5: Radiation Increments
       SW_INCS, LW_INCS, &
! PAR radiation increment
       DIRPAR, &
      ! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP, &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3, &
!  tropopause-based ozone
       TPPSOZONE, &
      ! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK, &       
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6, &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3, &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD, &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, &
       CO2, CH4_STOCH, O3_STOCH, &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5, &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10, &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15, &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20, &
! CABLE
!       TSOIL_TILE, SMCL_TILE, STHU_TILE, STHF_TILE, SNOW_DEPTH3L,       &
       TSOIL_TILE, SMCL_TILE,  STHF_TILE, SNOW_DEPTH3L,                 &
       SNOW_MASS3L, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SNOW_AGE,       & 
       SNOW_FLG3L,                                                      &
! Lestevens Sept 2012 - adding progs for CASACNP
       CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
       NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
       cable_lai,                                                       &

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC, &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC, &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC, &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC, TRACER_LBC, &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND, &
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND, &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND, &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND, &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND, &
       MURK_LBC_TEND, TRACER_LBC_TEND, &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB, &
      ! 2.2: Data variables stored in secondary space.
       PSTAR, &
      ! 2.3: Cloud fields
       CCB, CCT, CCLWP, &
      ! 2.4: Boundary layer fields
       ZH, & 
      ! Standard deviation of turbulent fluctuations of layer 1
       T1_SD, &
      ! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD, &
      ! Number of model levels in the  turbulently mixed layer
       NTML, &
      ! Top level for turb mixing in any decoupled Sc layer
       NTDSC, &
      ! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS, & 
      ! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT, &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN, &
      ! 2.5: Vegetation Ancillary fields
       CANOPY_WATER, SURF_CAP, SURF_RESIST, ROOT_DEPTH, INFILT, &
       VEG_FRAC, LAI, CANHT, Z0, SFA, MDSA, GS, &
      ! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2, &
       OROG_GRAD_X, OROG_GRAD_Y, &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY, &
      ! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS, &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, &
      ! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT, &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND, &
       DUST_MREL1, DUST_MREL2, DUST_MREL3, &
       DUST_MREL4, DUST_MREL5, DUST_MREL6, &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM, &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX, &
      ! Tracer Fluxes - kdcorbin, 05/10
       TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4, &
       TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8, &
       TRACER_FLUX9, TRACER_FLUX10, TRACER_FLUX11, TRACER_FLUX12, &
       TRACER_FLUX13, TRACER_FLUX14, TRACER_FLUX15, TRACER_FLUX16, &
       TRACER_FLUX17, TRACER_FLUX18, TRACER_FLUX19, TRACER_FLUX20, &

      ! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5, &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10, &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15, &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20, &
      !   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX, &
      !   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5, &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9, &
       LAI_PFT, CANHT_PFT, DISTURB_VEG, &
       SOIL_ALB, SOIL_CARB, &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4, &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC, &
       RSP_W_PFT_ACC, RSP_S_ACC, &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4, &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE, &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, &
       DOLR_FIELD, &
       LW_DOWN, SW_TILE_RTS, &
!! MODIFIED BY AT
!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!       TSLAB, TCLIM, HCLIM, CHEAT, OIFLX, UICE, VICE, &
!       SIG11NE, SIG11SE, SIG11SW, SIG11NW, &
!       SIG12NE, SIG12SE, SIG12SW, SIG12NW, &
!       SIG22NE, SIG22SE, SIG22SW, SIG22NW, &
!! END MODIFIED BY AT
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS, &
!   2.15: Fields carried forward from previous version.
!         May not be required
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho, &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP, &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND, &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET, &
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE, &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM, &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT, &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN, &
       C_SOLAR,C_BLUE,C_DOWN,C_LONGWAVE,C_TAUX,C_TAUY,C_WINDMIX, &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_BOTMELTN,C_TOPMELTN,C_LSRAIN, &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_PRESS, C_U10, C_V10, &
! UKCA oxidant fields
        OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA, & 
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR, &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK, &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4, &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC, &
       ARCLDLTA_DL, &
! Convective Cloud Fields
       LCBASE, CCW_RAD, &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM, &

! End arg_atm_fields.h
! ARGCONA start
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field sin_u_latitude. J F Thomson.
!  5.0   02/06/99  Insert C-P C-grid constants. M L Gallani.
!  5.3   01/10/01  Add fields for chequerboard radiation. S Cusack
! 6.1  04/08/04  Add separate arrays diff_coeff_u, diff_coeff_v
!                                                     Terry Davies
! 6.2  05/01/06   Add true_latitude. Yongming Tang
! 6.2  14/12/06  Add separate arrays VarRes Array co-ordinates
!                                                       Terry Davies

! argcona.h contained constants for the atmosphere.
! As of vn6.6 these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, AD_MASK_TROP_MOD, ROT_COEFF_MOD

! ARGLNDM Constants for physics routines
     &  land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     & ErrorStatus,CMessage)

      ENDIF ! Diagnostics required for this section

      IF (L_pc2 .and. ErrorStatus == 0 .and. l_physics) then

        deallocate(t_inc_pres)
        deallocate(q_inc_pres)
        deallocate(qcl_inc_pres)
        deallocate(qcf_inc_pres)
        deallocate(cf_inc_pres)
        deallocate(cfl_inc_pres)
        deallocate(cff_inc_pres)
        deallocate(t_dini)
        deallocate(q_dini)
        deallocate(qcl_dini)
        deallocate(qcf_dini)
        deallocate(cf_dini)
        deallocate(cfl_dini)
        deallocate(cff_dini)

      End If  ! L_pc2 and ErrorStatus == 0

! section 30: climate diagnostics
      IF(      SF(0,30)                                                 &
                                ! Diagnostics required for this section
     &   .AND. ErrorStatus == 0) THEN
! size of diagnostic space
        ALLOCATE (STASHwork30(STASH_maxlen(30,A_im)))
! DEPENDS ON: st_diag3
        CALL St_diag3(STASHwork30,STASH_maxlen(30,A_im),                &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
     &  A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
     &  A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
     &  A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Start arg_atm_fields.h

! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!
! Current Code Owner: A. Treshansky
!


      ! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP, &
      ! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV, &
      ! 1.2: Data variables stored in secondary space.
       P, & 
      ! Pressure on theta levels
       P_THETA_LEVELS, &
      ! Exner pressure on theta levels
       EXNER_THETA_LEVELS, &
      ! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN, &
      ! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF, &
      ! 1.5: Radiation Increments
       SW_INCS, LW_INCS, &
! PAR radiation increment
       DIRPAR, &
      ! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP, &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3, &
!  tropopause-based ozone
       TPPSOZONE, &
      ! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK, &       
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6, &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3, &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD, &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, &
       CO2, CH4_STOCH, O3_STOCH, &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5, &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10, &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15, &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20, &
! CABLE
!       TSOIL_TILE, SMCL_TILE, STHU_TILE, STHF_TILE, SNOW_DEPTH3L,       &
       TSOIL_TILE, SMCL_TILE,  STHF_TILE, SNOW_DEPTH3L,                 &
       SNOW_MASS3L, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SNOW_AGE,       & 
       SNOW_FLG3L,                                                      &
! Lestevens Sept 2012 - adding progs for CASACNP
       CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
       NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
       cable_lai,                                                       &

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC, &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC, &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC, &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC, TRACER_LBC, &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND, &
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND, &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND, &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND, &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND, &
       MURK_LBC_TEND, TRACER_LBC_TEND, &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB, &
      ! 2.2: Data variables stored in secondary space.
       PSTAR, &
      ! 2.3: Cloud fields
       CCB, CCT, CCLWP, &
      ! 2.4: Boundary layer fields
       ZH, & 
      ! Standard deviation of turbulent fluctuations of layer 1
       T1_SD, &
      ! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD, &
      ! Number of model levels in the  turbulently mixed layer
       NTML, &
      ! Top level for turb mixing in any decoupled Sc layer
       NTDSC, &
      ! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS, & 
      ! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT, &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN, &
      ! 2.5: Vegetation Ancillary fields
       CANOPY_WATER, SURF_CAP, SURF_RESIST, ROOT_DEPTH, INFILT, &
       VEG_FRAC, LAI, CANHT, Z0, SFA, MDSA, GS, &
      ! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2, &
       OROG_GRAD_X, OROG_GRAD_Y, &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY, &
      ! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS, &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, &
      ! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT, &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND, &
       DUST_MREL1, DUST_MREL2, DUST_MREL3, &
       DUST_MREL4, DUST_MREL5, DUST_MREL6, &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM, &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX, &
      ! Tracer Fluxes - kdcorbin, 05/10
       TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4, &
       TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8, &
       TRACER_FLUX9, TRACER_FLUX10, TRACER_FLUX11, TRACER_FLUX12, &
       TRACER_FLUX13, TRACER_FLUX14, TRACER_FLUX15, TRACER_FLUX16, &
       TRACER_FLUX17, TRACER_FLUX18, TRACER_FLUX19, TRACER_FLUX20, &

      ! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5, &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10, &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15, &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20, &
      !   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX, &
      !   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5, &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9, &
       LAI_PFT, CANHT_PFT, DISTURB_VEG, &
       SOIL_ALB, SOIL_CARB, &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4, &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC, &
       RSP_W_PFT_ACC, RSP_S_ACC, &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4, &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE, &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, &
       DOLR_FIELD, &
       LW_DOWN, SW_TILE_RTS, &
!! MODIFIED BY AT
!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!       TSLAB, TCLIM, HCLIM, CHEAT, OIFLX, UICE, VICE, &
!       SIG11NE, SIG11SE, SIG11SW, SIG11NW, &
!       SIG12NE, SIG12SE, SIG12SW, SIG12NW, &
!       SIG22NE, SIG22SE, SIG22SW, SIG22NW, &
!! END MODIFIED BY AT
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS, &
!   2.15: Fields carried forward from previous version.
!         May not be required
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho, &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP, &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND, &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET, &
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE, &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM, &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT, &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN, &
       C_SOLAR,C_BLUE,C_DOWN,C_LONGWAVE,C_TAUX,C_TAUY,C_WINDMIX, &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_BOTMELTN,C_TOPMELTN,C_LSRAIN, &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_PRESS, C_U10, C_V10, &
! UKCA oxidant fields
        OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA, & 
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR, &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK, &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4, &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC, &
       ARCLDLTA_DL, &
! Convective Cloud Fields
       LCBASE, CCW_RAD, &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM, &

! End arg_atm_fields.h
! ARGCONA start
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field sin_u_latitude. J F Thomson.
!  5.0   02/06/99  Insert C-P C-grid constants. M L Gallani.
!  5.3   01/10/01  Add fields for chequerboard radiation. S Cusack
! 6.1  04/08/04  Add separate arrays diff_coeff_u, diff_coeff_v
!                                                     Terry Davies
! 6.2  05/01/06   Add true_latitude. Yongming Tang
! 6.2  14/12/06  Add separate arrays VarRes Array co-ordinates
!                                                       Terry Davies

! argcona.h contained constants for the atmosphere.
! As of vn6.6 these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, AD_MASK_TROP_MOD, ROT_COEFF_MOD

! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &    energy_corr_now,                                              &
     &    inc_u, inc_v, inc_w, inc_t,                                   &
     &    inc_q, inc_qcl, inc_qcf,                                      &
!PC2     &    inc_cf, inc_cfl, inc_cff,
     &    inc_rho,                                                      &
     &    ErrorStatus,CMessage)

        DEALLOCATE (STASHwork30) ! Clear space
        DEALLOCATE (inc_q)
        DEALLOCATE (inc_qcl)
        DEALLOCATE (inc_qcf)
        DEALLOCATE (inc_cf)
        DEALLOCATE (inc_cfl)
        DEALLOCATE (inc_cff)
      ENDIF ! Diagnostics required for this section

      If( ErrorStatus == 0) then
        DEALLOCATE (inc_rho)
        DEALLOCATE (inc_t)
        DEALLOCATE (inc_u)
        DEALLOCATE (inc_v)
        DEALLOCATE (inc_w)
      Endif

! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport 
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

! DEPENDS ON: timer
      IF(LTIMER) CALL TIMER('AS9 End TStep Diags',6)

! section 0: extraction of primary variables

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,0,STASHWORK0_dummy,                          &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &           ErrorStatus,CMessage)
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',4)
! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,33,STASHWORK0_dummy,                         &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &           ErrorStatus,CMessage)
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',4)
! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF
! Section 34: extraction of UKCA tracer variables
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,34,STASHWORK0_dummy,                         &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &           ErrorStatus,CMessage)
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',4)

! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

! ---------------------------------------------------------------
! Section 9.1  Optional diagnostic printing of w, divergence
!               lapse rate and bottom level theta
! ---------------------------------------------------------------
       if( L_diag_print .and. L_diag_print_ops ) then
! DEPENDS ON: print_ops_diag
         Call Print_ops_diag(                                           &
     &                  W, THETA, row_length, rows,                     &
     &                  model_levels, model_domain,                     &
     &                  global_row_length, global_rows,                 &
     &                  offx, offy, halo_i, halo_j,                     &
     &                  mype, nproc, gc_proc_row_group,                 &
     &                  timestep_number, print_step, diag_interval,     &
     &                  rpemax, rpemin, ipesum,                         &
     &                  L_print_pe, L_print_wmax, L_print_theta1,       &
     &                  max_w_run, min_theta1_run,                      &
     &                  time_w_max, time_theta1_min )
! DEPENDS ON: um_fort_flush
         if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
       elseif( L_diag_print ) then
! DEPENDS ON: print_diag
         Call Print_diag(                                               &
     &                  U, V, THETA,                                    &
     &                  RHO, W,                                         &
     &                  Q, QCL, QCF,                                    &
     &                  rows, n_rows, row_length,                       &
     &                  model_levels, wet_levels, model_domain,         &
     &                  global_row_length, global_rows,                 &
     &                  r_theta_levels, r_rho_levels,                   &
     &                  r_at_u, r_at_v , Pi,                            &
     &                  FV_sec_theta_latitude, FV_cos_theta_latitude,   &
     &                  cos_v_latitude,                                 &
     &                  cos_theta_longitude, sin_theta_longitude,       &
     &                  offx, offy, halo_i, halo_j,                     &
     &                  mype, nproc, at_extremity, datastart,           &
     &                  gc_proc_row_group, delta_lambda, delta_phi,     &
     &                  timestep_number, print_step, diag_interval,     &
     &                  rpemax, rpemin, ipesum, rpesum, w_print_limit,  &
     &                  L_print_pe, L_print_w,                          &
     &                  L_print_wmax, L_print_max_wind,                 &
     &                  L_print_div, L_print_lapse, L_print_theta1,     &
     &                  L_print_shear, L_diag_wind, L_diag_noise,       &
     &                  max_w_run, max_wind_run, min_theta1_run,        &
     &                  dtheta1_run, max_div_run, min_div_run,          &
     &                  min_lapse_run, max_shear_run, time_max_shear,   &
     &                  time_div_max, time_div_min, time_lapse_min,     &
     &                  time_w_max, time_max_wind, time_theta1_min,     &
     &                  max_KE_run, min_KE_run, max_noise_run,          &
     &                  time_KE_max, time_KE_min, time_noise_max )
! DEPENDS ON: um_fort_flush
         if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
       endif     !  L_diag_print

! Logical first_atmstep_call is true on the first call to ATM_STEP
! and is set to false at the end of ATM_STEP (uses the SAVE command)

        first_atmstep_call = .false.

      RETURN
      END SUBROUTINE Atm_Step

