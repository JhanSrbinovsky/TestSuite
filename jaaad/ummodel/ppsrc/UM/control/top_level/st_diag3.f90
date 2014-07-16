
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ End of Timestep climate diagnostics control routine
!
! Subroutine Interface:
      SUBROUTINE ST_DIAG3( STASHWORK,INT30,                             &
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
     &    inc_q, inc_qcl, inc_qcf, inc_rho,                             &
     &                    ICODE,CMESSAGE)

      Use level_heights_mod
      Use trignometric_mod,  Only : sin_theta_latitude,                 &
     &     cos_theta_latitude, sin_theta_longitude, cos_theta_longitude
      Use rad_mask_trop_mod, Only : max_trop_level,min_trop_level

      IMPLICIT NONE
! Description:
! Calculates end of timestep diagnostics (held in STASH section 30)
! of which there are three types:
! 1) Model level fields and products
! 2) Pressure level fields and products ie U V T UU VV TT UT VT
! 3) Single level fields eg mountain torque and column integrations
!    eq atmospheric mass
!
! Method:
!   This routine sets controls the EoT diagnostics, and sets up
!   array sizes etc for a call to EOT_DIAG which calculates the
!   fields. STASHWORK is filled with the required data, and
!   passed up to atm_step2, ready for a call to stash.
!   The stash items are split into centuries, as follows:
!
!   000s   Standard Model level fields and products
!   100s   Other Model level fields
!   200s   Standard Pressure level fields and products
!   300s   Other Pressure level fields
!   400s   Surface fields and Column integrals
!
!  Standard diagnostics are assigned a number, which are:
!  1 U wind, 2 V wind, 3 W wind, 4 Temperature, 5 Q, 6 RH, 7 Z, 8 Omega
!
!
!  Further user and technical documentation can be found under:
!  http://www-hc/~hadsu/
!
! Current Code Owner: Simon Wilson
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.1     03/03/2000 Original code. SS Wilson
!
! 5.2     15/08/2000 Pstar added. SS Wilson
!
! 5.2     22/09/2000 Extra diagnostics and bug fixes. SS Wilson
!
! 5.3     08/08/2001 Tropopause diagnostics added. DGH Tan
! 5.4     24/04/2002 Add RH on model levels, wbig and the square of
!                    increments. R A Stratton
!         04/07/2002 Add column integral of cvT and gr. R A Stratton.
! 5.5     13/12/2002 Add wbig=0.1, Correct an extend model level
!                    products. R A Stratton
!         18/12/2002 Add more column integrals. R A Stratton.
!         31/12/2002 Dry mass weighting on model levels. R A Stratton
! 5.5     23/01/2003 Correct mountain torque diagnostic and add
!                    surface pressure drag diagnostics. S. Webster
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
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
! C_ETA_PMSL start
! ETA_PMSL is the ETA value used to determine the model level
! used in PMSL reduction calculations.
      REAL,PARAMETER:: ETA_PMSL=0.795
      INTEGER LEVNO_PMSL_CALC   !  Model level for PMSL reduction calc.
      COMMON /PMSLCALC/ LEVNO_PMSL_CALC
! C_ETA_PMSL end

      Integer, parameter ::                                             &
     &  npress_diags=8                                                  &
                                 ! No. of diags on pressure levels
     &, nmodel_diags=7           ! No. of diags on modellevels

!   Scalar  arguments with intent(in):

!   Array  arguments with intent(in):
      real energy_corr_now                                              &
     &,  inc_u(1-offx:row_length+offx, 1-offy:rows+offy,                &
     &        model_levels)                                             &
     &, inc_v(1-offx:row_length+offx, 1-offy:n_rows+offy,               &
     &        model_levels)                                             &
     &, inc_w(row_length, rows, model_levels)                           &
     &, inc_t(1-offx:row_length+offx,                                   &
     &               1-offy:rows+offy, model_levels)                    &
     &, inc_q(1-halo_i:row_length+halo_i,                               &
     &               1-halo_j:rows+halo_j,wet_levels)                   &
     &, inc_qcl(1-halo_i:row_length+halo_i,                             &
     &               1-halo_j:rows+halo_j,wet_levels)                   &
     &, inc_qcf(1-halo_i:row_length+halo_i,                             &
     &               1-halo_j:rows+halo_j,wet_levels)                   &
     &, inc_rho(1-offx:row_length+offx,                                 &
     &             1-offy:rows+offy, model_levels)

!   Scalar arguments with intent(InOut):
      INTEGER                                                           &
     &        INT30,                                                    &
                                ! Dummy for STASH_MAXLEN(30)
     &        ICODE              ! Out return code : 0 Normal exit
                                !                 : >0 Error exit

      CHARACTER*256                                                     &
     &        CMESSAGE          ! Out error message if ICODE > 0

!   Array  arguments with intent(InOut):
      REAL                                                              &
     &  STASHWORK(*)

!   Scalar arguments with intent(out):
      INTEGER                                                           &
     &  U_M_LEVS,                                                       &
     &  V_M_LEVS,                                                       &
     &  W_M_LEVS,                                                       &
     &  T_M_LEVS,                                                       &
     &  Q_M_LEVS,                                                       &
     &  Z_M_LEVS,                                                       &
     &  KE_M_LEVS,                                                      &
     &  wbig_m_levs,                                                    &
     &  wbig2_m_levs,                                                   &
     &  dry_mass_m_levs,                                                &
     &  Rh_m_levs,                                                      &
     &  U_MM_LEVS,                                                      &
     &  V_MM_LEVS,                                                      &
     &  W_MM_LEVS,                                                      &
     &  T_MM_LEVS,                                                      &
     &  Q_MM_LEVS,                                                      &
     &  Z_MM_LEVS,                                                      &
     &  KE_MM_LEVS,                                                     &
     &  TEOT_M_LEVS,                                                    &
     &  U_INC_LEVS,                                                     &
     &  V_INC_LEVS,                                                     &
     &  W_INC_LEVS,                                                     &
     &  T_INC_LEVS,                                                     &
     &  Q_INC_LEVS,                                                     &
     &  QCL_INC_LEVS,                                                   &
     &  QCF_INC_LEVS,                                                   &
     &  RHO_INC_LEVS,                                                   &
     &  U2_INC_LEVS,                                                    &
     &  V2_INC_LEVS,                                                    &
     &  W2_INC_LEVS,                                                    &
     &  T2_INC_LEVS,                                                    &
     &  Q2_INC_LEVS,                                                    &
     &  QCL2_INC_LEVS,                                                  &
     &  QCF2_INC_LEVS,                                                  &
     &  RHO2_INC_LEVS,                                                  &
     &  U_P_LEVS,                                                       &
     &  V_P_LEVS,                                                       &
     &  W_P_LEVS,                                                       &
     &  T_P_LEVS,                                                       &
     &  Q_P_LEVS,                                                       &
     &  RH_P_LEVS,                                                      &
     &  Z_P_LEVS,                                                       &
     &  OM_P_LEVS,                                                      &
     &  HEAVY_P_LEVS,                                                   &
     &  TV_P_LEVS,                                                      &
     &  TVOM_P_LEVS,                                                    &
     &  FIELD1_P_LEVS,                                                  &
     &  FIELD2_P_LEVS,                                                  &
     &  n_levels                                                        &
     &  ,im_ident                                                       &
                                !  Internal Model Identifier
     &  ,im_index               !  Internal Model Index for Stash Arrays

      INTEGER                                                           &
     &  PT001,PT002,PT003,PT004,PT005,PT006,PT007,                      &
     &  PT101,PT102,PT103,PT104,PT105,PT106,PT107,                      &
     &  PT111,PT112,PT113,PT114,PT115,                                  &
     &  PT181,PT182,PT183,PT184,PT185,PT186,PT187,PT188,                &
     &  PT171,PT172,PT173,PT174,PT175,PT176,PT177,PT178,                &
     &  PT201,PT202,PT203,PT204,PT205,PT206,PT207,PT208,                &
     &  PT301,PT302,PT303,                                              &
     &  PT401,PT402,PT403,PT404,PT405,PT406,PT407,PT408,PT409,PT410,    &
     &  PT411,PT412,PT413,PT414,PT415,PT416,PT417,PT418,PT419,PT420,    &
     &  PT421,PT422,PT423,PT424,PT425,PT426,PT427,PT428,PT429,PT430,    &
     &  PT431,PT432,PT433,PT434,PT435,PT436,PT437,PT438,PT439           &
     & ,PT440,PT441,PT451,PT452,PT453,PT454

!   Array  arguments with intent(out):
      REAL                                                              &
     &  U_PRESS(NUM_STASH_LEVELS)                                       &
     &  ,V_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,W_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,T_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,Q_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,RH_PRESS(NUM_STASH_LEVELS)                                     &
     &  ,Z_PRESS(NUM_STASH_LEVELS)                                      &
                                   ! Z pressure levels
     &  ,OM_PRESS(NUM_STASH_LEVELS)                                     &
     &  ,HEAVY_PRESS(NUM_STASH_LEVELS)                                  &
     &  ,PRESS_LEVS(NUM_STASH_LEVELS)

      LOGICAL U_M_LIST(MODEL_LEVELS)                                    &
     &  ,V_M_LIST(MODEL_LEVELS)                                         &
     &  ,W_M_LIST(MODEL_LEVELS)                                         &
     &  ,T_M_LIST(MODEL_LEVELS)                                         &
     &  ,Q_M_LIST(MODEL_LEVELS)                                         &
     &  ,Z_M_LIST(MODEL_LEVELS)                                         &
     &  ,KE_M_LIST(MODEL_LEVELS)                                        &
     &  ,wbig_m_list(MODEL_LEVELS)                                      &
     &  ,wbig2_m_list(MODEL_LEVELS)                                     &
     &  ,RH_m_list(MODEL_LEVELS)                                        &
     &  ,DRY_MASS_M_LIST(MODEL_LEVELS)                                  &
     &  ,U_MM_LIST(MODEL_LEVELS)                                        &
     &  ,V_MM_LIST(MODEL_LEVELS)                                        &
     &  ,W_MM_LIST(MODEL_LEVELS)                                        &
     &  ,T_MM_LIST(MODEL_LEVELS)                                        &
     &  ,Q_MM_LIST(MODEL_LEVELS)                                        &
     &  ,Z_MM_LIST(MODEL_LEVELS)                                        &
     &  ,KE_MM_LIST(MODEL_LEVELS)                                       &
     &  ,TEOT_M_LIST(MODEL_LEVELS)                                      &
     &  ,U_INC_LIST(MODEL_LEVELS)                                       &
     &  ,V_INC_LIST(MODEL_LEVELS)                                       &
     &  ,W_INC_LIST(MODEL_LEVELS)                                       &
     &  ,T_INC_LIST(MODEL_LEVELS)                                       &
     &  ,Q_INC_LIST(MODEL_LEVELS)                                       &
     &  ,QCL_INC_LIST(MODEL_LEVELS)                                     &
     &  ,QCF_INC_LIST(MODEL_LEVELS)                                     &
     &  ,RHO_INC_LIST(MODEL_LEVELS)                                     &
     &  ,U2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,V2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,W2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,T2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,Q2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,QCL2_INC_LIST(MODEL_LEVELS)                                    &
     &  ,QCF2_INC_LIST(MODEL_LEVELS)                                    &
     &  ,RHO2_INC_LIST(MODEL_LEVELS)

      logical prod_m_list(model_levels,nmodel_diags,nmodel_diags)

      INTEGER                                                           &
     &  prod_IND(NUM_STASH_LEVELS*2,npress_diags,npress_diags),         &
     &  prod_p_levs(npress_diags,npress_diags)

! Local parameters:

! Local scalars:
      INTEGER                                                           &
     &  I,                                                              &
     &  NI,                                                             &
     &  K,                                                              &
     &  ISL,                                                            &
     &  BL,                                                             &
     &  TL,                                                             &
     &  LEVEL,                                                          &
     &  code30,nii,nij,j

      character*2 cfield1,cfield2

! Local dynamic arrays:
      REAL                                                              &
     &  FIELD1_PRESS(NUM_STASH_LEVELS)                                  &
     &  ,FIELD2_PRESS(NUM_STASH_LEVELS)

! Function & Subroutine calls:

      EXTERNAL                                                          &
     &  STASH,                                                          &
     &  TIMER,                                                          &
     &  EOT_DIAG,                                                       &
     &  SET_LEVELS_LIST,                                                &
     &  CHECK_PROD_LEVS

!- End of header

!     Set to atmosphere internal model
      im_ident = atmos_im
      im_index = internal_model_index(im_ident)



!--------------------Extract Reqd Model Levels for U-------------
      ISL=STINDEX(1,001,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U_M_LIST(I)) U_M_LEVS=U_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          U_M_LIST(I)=.FALSE.
        ENDDO
        U_M_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for V-------------
      ISL=STINDEX(1,002,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V_M_LIST(I)) V_M_LEVS=V_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          V_M_LIST(I)=.FALSE.
        ENDDO
        V_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for W-------------
      ISL=STINDEX(1,003,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W_M_LIST(I)) W_M_LEVS=W_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          W_M_LIST(I)=.FALSE.
        ENDDO
        W_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T-------------
      ISL=STINDEX(1,004,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T_M_LIST(I)) T_M_LEVS=T_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          T_M_LIST(I)=.FALSE.
        ENDDO
        T_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q-------------
      ISL=STINDEX(1,005,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q_M_LIST(I)) Q_M_LEVS=Q_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Q_M_LIST(I)=.FALSE.
        ENDDO
        Q_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Z-------------
      ISL=STINDEX(1,006,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Z_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Z_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Z_M_LIST(I)) Z_M_LEVS=Z_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Z_M_LIST(I)=.FALSE.
        ENDDO
        Z_M_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for KE-------------
      ISL=STINDEX(1,007,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    KE_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        KE_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (KE_M_LIST(I)) KE_M_LEVS=KE_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          KE_M_LIST(I)=.FALSE.
        ENDDO
        KE_M_LEVS=1
      ENDIF
!--------------Extract Reqd Model Levels for Dry mass weighting----
      ISL=STINDEX(1,115,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &  dry_mass_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        DRY_MASS_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (DRY_MASS_M_LIST(I)) DRY_MASS_M_LEVS=DRY_MASS_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          DRY_MASS_M_LIST(I)=.FALSE.
        ENDDO
        DRY_MASS_M_LEVS=1
      ENDIF

!-Extract Reqd model levels for products and store in prod_m_list

      do i=1,nmodel_diags
        do j=1,nmodel_diags
          code30=i*10+j
          IF (SF(code30,30)) THEN
            ISL=STINDEX(1,code30,30,im_index)
! DEPENDS ON: set_levels_list
            CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),     &
     &        PROD_M_LIST(1,i,j),                                       &
     &        STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
          ELSE
            DO k=1,MODEL_LEVELS
              prod_M_LIST(k,i,j)=.FALSE.
            ENDDO
          ENDIF
        enddo
      enddo


!--------------------Extract Reqd Model Levels for wbig----------
      ISL=STINDEX(1,112,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Wbig_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Wbig_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Wbig_M_LIST(I)) Wbig_M_LEVS=Wbig_M_LEVS+1
        ENDDO
      ELSE
        Wbig_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for wbig2----------
      ISL=STINDEX(1,114,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Wbig2_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Wbig2_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Wbig2_M_LIST(I)) Wbig2_M_LEVS=Wbig2_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Wbig2_M_LIST(I)=.FALSE.
        ENDDO
         Wbig2_M_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for RH ----------
      ISL=STINDEX(1,113,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    RH_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        RH_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (RH_M_LIST(I)) RH_M_LEVS=RH_M_LEVS+1
        ENDDO
      ELSE
        RH_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for U mass weighted
      ISL=STINDEX(1,101,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U_MM_LIST(I)) U_MM_LEVS=U_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          U_MM_LIST(I)=.FALSE.
        ENDDO
        U_MM_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for V mass weighted
      ISL=STINDEX(1,102,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V_MM_LIST(I)) V_MM_LEVS=V_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          V_MM_LIST(I)=.FALSE.
        ENDDO
        V_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for W mass weighted
      ISL=STINDEX(1,103,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W_MM_LIST(I)) W_MM_LEVS=W_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          W_MM_LIST(I)=.FALSE.
        ENDDO
        W_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T mass weighted
      ISL=STINDEX(1,104,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T_MM_LIST(I)) T_MM_LEVS=T_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          T_MM_LIST(I)=.FALSE.
        ENDDO
        T_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q mass weighted
      ISL=STINDEX(1,105,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q_MM_LIST(I)) Q_MM_LEVS=Q_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Q_MM_LIST(I)=.FALSE.
        ENDDO
        Q_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Z mass weighted
      ISL=STINDEX(1,106,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Z_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Z_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Z_MM_LIST(I)) Z_MM_LEVS=Z_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Z_MM_LIST(I)=.FALSE.
        ENDDO
        Z_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for KE mass weighted
      ISL=STINDEX(1,107,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    KE_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        KE_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (KE_MM_LIST(I)) KE_MM_LEVS=KE_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          KE_MM_LIST(I)=.FALSE.
        ENDDO
        KE_MM_LEVS=1
      ENDIF


!--------------------Extract Reqd Model Levels for T at EOT
      ISL=STINDEX(1,111,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    TEOT_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        TEOT_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (TEOT_M_LIST(I)) TEOT_M_LEVS=TEOT_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          TEOT_M_LIST(I)=.FALSE.
        ENDDO
        TEOT_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T increment
      ISL=STINDEX(1,181,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T_INC_LIST(I)) T_INC_LEVS=T_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          T_INC_LIST(I)=.FALSE.
        ENDDO
        T_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q increment
      ISL=STINDEX(1,182,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q_INC_LIST(I)) Q_INC_LEVS=Q_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Q_INC_LIST(I)=.FALSE.
        ENDDO
        Q_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCL increment
      ISL=STINDEX(1,183,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCL_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCL_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCL_INC_LIST(I)) QCL_INC_LEVS=QCL_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          QCL_INC_LIST(I)=.FALSE.
        ENDDO
        QCL_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCF increment
      ISL=STINDEX(1,184,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCF_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCF_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCF_INC_LIST(I)) QCF_INC_LEVS=QCF_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          QCF_INC_LIST(I)=.FALSE.
        ENDDO
        QCF_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for U increment
      ISL=STINDEX(1,185,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U_INC_LIST(I)) U_INC_LEVS=U_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          U_INC_LIST(I)=.FALSE.
        ENDDO
        U_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for V increment
      ISL=STINDEX(1,186,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V_INC_LIST(I)) V_INC_LEVS=V_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          V_INC_LIST(I)=.FALSE.
        ENDDO
        V_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for W increment
      ISL=STINDEX(1,187,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W_INC_LIST(I)) W_INC_LEVS=W_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          W_INC_LIST(I)=.FALSE.
        ENDDO
        W_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for RHO increment
      ISL=STINDEX(1,188,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    RHO_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        RHO_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (RHO_INC_LIST(I)) RHO_INC_LEVS=RHO_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          RHO_INC_LIST(I)=.FALSE.
        ENDDO
        RHO_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T increment **2
      ISL=STINDEX(1,171,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T2_INC_LIST(I)) T2_INC_LEVS=T2_INC_LEVS+1
        ENDDO
      ELSE
        T2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q increment **2
      ISL=STINDEX(1,172,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q2_INC_LIST(I)) Q2_INC_LEVS=Q2_INC_LEVS+1
        ENDDO
      ELSE
        Q2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCL increment **2
      ISL=STINDEX(1,173,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCL2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCL2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCL2_INC_LIST(I)) QCL2_INC_LEVS=QCL2_INC_LEVS+1
        ENDDO
      ELSE
        QCL2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCF increment **2
      ISL=STINDEX(1,174,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCF2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCF2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCF2_INC_LIST(I)) QCF2_INC_LEVS=QCF2_INC_LEVS+1
        ENDDO
      ELSE
        QCF2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for U increment **2
      ISL=STINDEX(1,175,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U2_INC_LIST(I)) U2_INC_LEVS=U2_INC_LEVS+1
        ENDDO
      ELSE
        U2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for V increment **2
      ISL=STINDEX(1,176,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V2_INC_LIST(I)) V2_INC_LEVS=V2_INC_LEVS+1
        ENDDO
      ELSE
        V2_INC_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for W increment **2
      ISL=STINDEX(1,177,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W2_INC_LIST(I)) W2_INC_LEVS=W2_INC_LEVS+1
        ENDDO
      ELSE
        W2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for RHO increment **2
      ISL=STINDEX(1,178,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    RHO2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        RHO2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (RHO2_INC_LIST(I)) RHO2_INC_LEVS=RHO2_INC_LEVS+1
        ENDDO
      ELSE
        RHO2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Pressures for U_COMP_P-------------

      ISL=STINDEX(1,201,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        U_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,U_P_LEVS
          U_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        U_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for V_COMP_P-------------

      ISL=STINDEX(1,202,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        V_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,V_P_LEVS
          V_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        V_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for W_COMP_P-------------

      ISL=STINDEX(1,203,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        W_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,W_P_LEVS
          W_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        W_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for T on wind grid-------

      ISL=STINDEX(1,204,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        T_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,T_P_LEVS
          T_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        T_P_LEVS=1
      END IF


!--------------------Extract Reqd Pressures for q on wind grid-------

      ISL=STINDEX(1,205,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        Q_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,Q_P_LEVS
          Q_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        Q_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for RH-----

      ISL=STINDEX(1,206,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        RH_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,RH_P_LEVS
          RH_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        RH_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for geopotential height-----

      ISL=STINDEX(1,207,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        Z_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,Z_P_LEVS
          Z_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        Z_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for Omega---

      ISL=STINDEX(1,208,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        OM_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,OM_P_LEVS
          OM_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        OM_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for products---


      do i=1,npress_diags
        do j=1,npress_diags
          if (i  ==  1) then
            cfield1='u'
            field1_press=u_press
            field1_p_levs=u_p_levs
          else if (i  ==  2) then
            cfield1='v'
            field1_press=v_press
            field1_p_levs=v_p_levs
          else if (i  ==  3) then
            cfield1='w'
            field1_press=w_press
            field1_p_levs=w_p_levs
          else if (i  ==  4) then
            cfield1='t'
            field1_press=t_press
            field1_p_levs=t_p_levs
          else if (i  ==  5) then
            cfield1='q'
            field1_press=q_press
            field1_p_levs=q_p_levs
          else if (i  ==  6) then
            cfield1='rh'
            field1_press=rh_press
            field1_p_levs=rh_p_levs
          else if (i  ==  7) then
            cfield1='z'
            field1_press=z_press
            field1_p_levs=z_p_levs
          else if (i  ==  8) then
            cfield1='om'
            field1_press=om_press
            field1_p_levs=om_p_levs
          endif

          if (j  ==  1) then
            cfield2='u'
            field2_press=u_press
            field2_p_levs=u_p_levs
          else if (j  ==  2) then
            cfield2='v'
            field2_press=v_press
            field2_p_levs=v_p_levs
          else if (j  ==  3) then
            cfield2='w'
            field2_press=w_press
            field2_p_levs=w_p_levs
          else if (j  ==  4) then
            cfield2='t'
            field2_press=t_press
            field2_p_levs=t_p_levs
          else if (j  ==  5) then
            cfield2='q'
            field2_press=q_press
            field2_p_levs=q_p_levs
          else if (j  ==  6) then
            cfield2='rh'
            field2_press=rh_press
            field2_p_levs=rh_p_levs
          else if (j  ==  7) then
            cfield2='z'
            field2_press=z_press
            field2_p_levs=z_p_levs
          else if (j  ==  8) then
            cfield2='om'
            field2_press=om_press
            field2_p_levs=om_p_levs
          endif
          code30=200+i*10+j
          IF (SF(code30,30)) THEN
            IF ((.NOT.SF(200+i,30)).OR.(.NOT.SF(200+j,30))) THEN
              CMESSAGE=' '
              CMESSAGE='ST_DIAG3 : '//cfield1//'*'//cfield2//           &
     &          ' error '//cfield1//' and '//cfield2//                  &
     &          ' must be requested'
              ICODE=1
              GOTO 999
            ELSE
              NI=-STLIST(10,STINDEX(1,code30,30,im_index))
              NII=-STLIST(10,STINDEX(1,200+i,30,im_index))
              NIJ=-STLIST(10,STINDEX(1,200+j,30,im_index))

! DEPENDS ON: check_prod_levs
              CALL CHECK_PROD_LEVS(                                     &
     &          NUM_STASH_LEVELS,STASH_LEVELS,NI,                       &
     &          FIELD1_P_LEVS,FIELD1_PRESS,                             &
     &          FIELD2_P_LEVS,FIELD2_PRESS,                             &
     &          prod_P_LEVS(i,j),prod_IND(1,i,j))
            END IF
          ELSE
            prod_P_LEVS(i,j)=1
          ENDIF
        enddo
      enddo

!     Do checks for other pressure fields

      if (sf(302,30).and..not.                                          &
     &  (sf(245,30).and.sf(204,30).and.sf(205,30))) then
        CMESSAGE='ST_DIAG3 : Virtual temperature'//                     &
     &          ' error T*Q must also be requested'
        ICODE=1
        GOTO 999
      endif

      if (sf(303,30).and..not.                                          &
     &  (sf(245,30).and.sf(204,30).and.                                 &
     &  sf(205,30).and.sf(208,30))) then
        CMESSAGE='ST_DIAG3 : Virtual temperature*omega'//               &
     &          ' error T*Q and omega must also be requested'
        ICODE=1
        GOTO 999
      endif



!--------------------Extract Reqd Pressures for Heavyside function---

      ISL=STINDEX(1,301,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        HEAVY_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,HEAVY_P_LEVS
          HEAVY_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        HEAVY_P_LEVS=1
      END IF

!Lset up level information for tv
      if (sf(302,30)) then
        if (sf(204,30).and.sf(205,30)) then
          tv_p_levs=prod_p_levs(4,5)
        else
          tv_p_levs=1
          CMESSAGE='ST_DIAG3 : T and Q must be selected to calculate '//&
     &      'virtual temperature'
          ICODE=1
          GOTO 999
        endif
      else
        tv_p_levs=1
      endif

!Lset up level information for tv*om
      if (sf(303,30)) then
        if (sf(204,30).and.sf(205,30).and.sf(208,30)) then
          tvom_p_levs=prod_p_levs(4,5)
        else
          tvom_p_levs=1
          CMESSAGE='ST_DIAG3 : T,Q and Omega'//                         &
     &      'must be selected to calculate '//                          &
     &      'virtual temperature'
          ICODE=1
          GOTO 999
        endif
      else
        tvom_p_levs=1
      endif



!-------------------Set up Pointers for STASHWORK -------------------
      PT001=SI(001,30,im_index)
      PT002=SI(002,30,im_index)
      PT003=SI(003,30,im_index)
      PT004=SI(004,30,im_index)
      PT005=SI(005,30,im_index)
      PT006=SI(006,30,im_index)
      PT007=SI(007,30,im_index)
      PT101=SI(101,30,im_index)
      PT102=SI(102,30,im_index)
      PT103=SI(103,30,im_index)
      PT104=SI(104,30,im_index)
      PT105=SI(105,30,im_index)
      PT106=SI(106,30,im_index)
      PT107=SI(107,30,im_index)
      PT111=SI(111,30,im_index)
      PT112=SI(112,30,im_index)
      PT113=SI(113,30,im_index)
      PT114=SI(114,30,im_index)
      PT115=SI(115,30,im_index)
      PT171=SI(171,30,im_index)
      PT172=SI(172,30,im_index)
      PT173=SI(173,30,im_index)
      PT174=SI(174,30,im_index)
      PT175=SI(175,30,im_index)
      PT176=SI(176,30,im_index)
      PT177=SI(177,30,im_index)
      PT178=SI(178,30,im_index)
      PT181=SI(181,30,im_index)
      PT182=SI(182,30,im_index)
      PT183=SI(183,30,im_index)
      PT184=SI(184,30,im_index)
      PT185=SI(185,30,im_index)
      PT186=SI(186,30,im_index)
      PT187=SI(187,30,im_index)
      PT188=SI(188,30,im_index)
      PT201=SI(201,30,im_index)
      PT202=SI(202,30,im_index)
      PT203=SI(203,30,im_index)
      PT204=SI(204,30,im_index)
      PT205=SI(205,30,im_index)
      PT206=SI(206,30,im_index)
      PT207=SI(207,30,im_index)
      PT208=SI(208,30,im_index)
      PT301=SI(301,30,im_index)
      PT302=SI(302,30,im_index)
      PT303=SI(303,30,im_index)
      PT401=SI(401,30,im_index)
      PT402=SI(402,30,im_index)
      PT403=SI(403,30,im_index)
      PT404=SI(404,30,im_index)
      PT405=SI(405,30,im_index)
      PT406=SI(406,30,im_index)
      PT407=SI(407,30,im_index)
      PT408=SI(408,30,im_index)
      PT409=SI(409,30,im_index)
      PT410=SI(410,30,im_index)
      PT411=SI(411,30,im_index)
      PT412=SI(412,30,im_index)
      PT413=SI(413,30,im_index)
      PT414=SI(414,30,im_index)
      PT415=SI(415,30,im_index)
      PT416=SI(416,30,im_index)
      PT417=SI(417,30,im_index)
      PT418=SI(418,30,im_index)
      PT419=SI(419,30,im_index)
      PT420=SI(420,30,im_index)
      PT421=SI(421,30,im_index)
      PT422=SI(422,30,im_index)
      PT423=SI(423,30,im_index)
      PT424=SI(424,30,im_index)
      PT425=SI(425,30,im_index)
      PT426=SI(426,30,im_index)
      PT427=SI(427,30,im_index)
      PT428=SI(428,30,im_index)
      PT429=SI(429,30,im_index)
      PT430=SI(430,30,im_index)
      PT431=SI(431,30,im_index)
      PT432=SI(432,30,im_index)
      PT433=SI(433,30,im_index)
      PT434=SI(434,30,im_index)
      PT435=SI(435,30,im_index)
      PT436=SI(436,30,im_index)
      PT437=SI(437,30,im_index)
      PT438=SI(438,30,im_index)
      PT439=SI(439,30,im_index)
      PT440=SI(440,30,im_index)
      PT441=SI(441,30,im_index)
      PT451=SI(451,30,im_index)
      PT452=SI(452,30,im_index)
      PT453=SI(453,30,im_index)
      PT454=SI(454,30,im_index)

! Initialise STASHWORK array because EOT_DIAG does not initialise halos
!* DIR$ CACHE_BYPASS STASHWORK
      DO I=1,INT30
        STASHWORK(I)=0.
      ENDDO

      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EOT_DIAG',3)
      END IF

! DEPENDS ON: eot_diag
      CALL Eot_diag(                                                    &
! Primary data: in
     &  PSTAR,P,RHO,U,V,W                                               &
     &  ,THETA,Q,QCL,QCF                                                &
     &  ,p_theta_levels ,exner_rho_levels                               &
     &  ,exner_theta_levels,energy_corr_now                             &
     &  ,inc_u, inc_v, inc_w, inc_t ,inc_q, inc_qcl, inc_qcf, inc_rho   &
! Grid sizes and definition: in
     &  ,rows,n_rows,row_length,model_levels,wet_levels,bl_levels       &
     &  ,theta_field_size,u_field_size,v_field_size                     &
     &  ,eta_theta_levels,eta_rho_levels,r_theta_levels,r_rho_levels    &
     &  ,sin_theta_latitude,cos_theta_latitude                          &
     &  ,sin_theta_longitude,cos_theta_longitude                        &
     &  ,Model_domain ,delta_lambda ,delta_phi,global_row_length        &
! Grid coordinates: in
     &  ,A_REALHD(rh_deltaEW),A_REALHD(rh_deltaNS)                      &
     &  ,A_REALHD(rh_baselat),A_REALHD(rh_baselong)                     &
     &  ,A_REALHD(rh_rotlat),A_REALHD(rh_rotlong)                       &
! Pre-calculated grid associated arrays: in
     &  ,r_at_u,r_at_v                                                  &
! Time information: in
     &  ,forecast_hrs                                                   &
! Pressure levels for output arrays: in
     &  ,U_PRESS,V_PRESS,W_PRESS,T_PRESS,Q_PRESS,RH_PRESS               &
     &  ,Z_PRESS,OM_PRESS,HEAVY_PRESS                                   &
! Flags to request each diagnostic output field: in
     &  ,SF(001,30),SF(002,30),SF(003,30),SF(004,30),SF(005,30)         &
     &  ,SF(006,30),SF(007,30)                                          &
     &  ,SF(101,30),SF(102,30),SF(103,30),SF(104,30),SF(105,30)         &
     &  ,SF(106,30),SF(107,30)                                          &
     &  ,SF(111,30),SF(112,30),SF(114,30),SF(113,30),SF(115,30)         &
     &  ,SF(181,30),SF(182,30),SF(183,30),SF(184,30),SF(185,30)         &
     &  ,SF(186,30),SF(187,30),SF(188,30)                               &
     &  ,SF(171,30),SF(172,30),SF(173,30),SF(174,30),SF(175,30)         &
     &  ,SF(176,30),SF(177,30),SF(178,30)                               &
     &  ,SF(201,30),SF(202,30),SF(203,30)                               &
     &  ,SF(204,30),SF(205,30),SF(206,30),SF(207,30),SF(208,30)         &
     &  ,SF(301,30),SF(302,30),SF(303,30)                               &
     &  ,SF(401,30),SF(402,30),SF(403,30),SF(404,30),SF(405,30)         &
     &  ,SF(406,30),SF(407,30),SF(408,30),SF(409,30),SF(410,30)         &
     &  ,SF(411,30),SF(412,30),SF(413,30),SF(414,30),SF(415,30)         &
     &  ,SF(416,30),SF(417,30),SF(418,30),SF(419,30),SF(420,30)         &
     &  ,SF(421,30),SF(422,30),SF(423,30),SF(424,30),SF(425,30)         &
     &  ,SF(426,30),SF(427,30),SF(428,30),SF(429,30),SF(430,30)         &
     &  ,SF(431,30),SF(432,30),SF(433,30),SF(434,30),SF(435,30)         &
     &  ,SF(436,30),SF(437,30),SF(438,30),SF(439,30)                    &
     &  ,SF(440,30),SF(441,30)                                          &
     &  ,SF(451,30),SF(452,30),SF(453,30),SF(454,30)                    &
! Diagnostics lengths: in
     &  ,u_m_levs,v_m_levs,w_m_levs,t_m_levs,q_m_levs,z_m_levs          &
     &  ,ke_m_levs,u_mm_levs,v_mm_levs,w_mm_levs,t_mm_levs,q_mm_levs    &
     &  ,z_mm_levs,ke_mm_levs,teot_m_levs,wbig_m_levs,wbig2_m_levs      &
     &  ,RH_m_levs,dry_mass_m_levs                                      &
     &  ,t_inc_levs,q_inc_levs,qcl_inc_levs,qcf_inc_levs                &
     &  ,u_inc_levs,v_inc_levs,w_inc_levs,rho_inc_levs                  &
     &  ,t2_inc_levs,q2_inc_levs,qcl2_inc_levs,qcf2_inc_levs            &
     &  ,u2_inc_levs,v2_inc_levs,w2_inc_levs,rho2_inc_levs              &
     &  ,u_p_levs,v_p_levs,w_p_levs, t_p_levs,q_p_levs,rh_p_levs        &
     &  ,z_p_levs,om_p_levs,heavy_p_levs,tv_p_levs,tvom_p_levs          &
     &  ,prod_p_levs                                                    &
                     ! pressure levels required for each product
     &  ,npress_diags,nmodel_diags                                      &
! Tropopause index bounds
     &  ,min_trop_level,max_trop_level                                  &
! Model levels indicies: in
     &  ,u_m_list,v_m_list,w_m_list,t_m_list,q_m_list,z_m_list,ke_m_list&
     &  ,u_mm_list,v_mm_list,w_mm_list,t_mm_list,q_mm_list,z_mm_list    &
     &  ,ke_mm_list,teot_m_list,wbig_m_list,wbig2_m_list,RH_m_list      &
     &  ,dry_mass_m_list                                                &
     &  ,t_inc_list,q_inc_list,qcl_inc_list,qcf_inc_list                &
     &  ,u_inc_list,v_inc_list,w_inc_list,rho_inc_list                  &
     &  ,t2_inc_list,q2_inc_list,qcl2_inc_list,qcf2_inc_list            &
     &  ,u2_inc_list,v2_inc_list,w2_inc_list,rho2_inc_list              &
     &  ,prod_m_list                                                    &
                     ! index of model levels required for product
! Product levels: in
     &  ,prod_ind                                                       &
                  ! index of pressure levels required for product
! Diagnostic arrays: out
     &  ,SI(1,30,im_index) ,Stashwork                                   &
     &  ,Stashwork(pt001),Stashwork(pt002),Stashwork(pt003)             &
     &  ,Stashwork(pt004),Stashwork(pt005),Stashwork(pt006)             &
     &  ,Stashwork(pt007)                                               &
     &  ,Stashwork(pt101),Stashwork(pt102),Stashwork(pt103)             &
     &  ,Stashwork(pt104),Stashwork(pt105),Stashwork(pt106)             &
     &  ,Stashwork(pt107),Stashwork(pt111),Stashwork(pt112)             &
     &  ,Stashwork(pt114),Stashwork(pt113),Stashwork(pt115)             &
     &  ,Stashwork(pt181),Stashwork(pt182),Stashwork(pt183)             &
     &  ,Stashwork(pt184),Stashwork(pt185),Stashwork(pt186)             &
     &  ,Stashwork(pt187),Stashwork(pt188)                              &
     &  ,Stashwork(pt171),Stashwork(pt172),Stashwork(pt173)             &
     &  ,Stashwork(pt174),Stashwork(pt175),Stashwork(pt176)             &
     &  ,Stashwork(pt177),Stashwork(pt178)                              &
     &  ,Stashwork(pt201),Stashwork(pt202),Stashwork(pt203)             &
     &  ,Stashwork(pt204),Stashwork(pt205),Stashwork(pt206)             &
     &  ,Stashwork(pt207),Stashwork(pt208)                              &
     &  ,Stashwork(pt301),Stashwork(pt302),Stashwork(pt303)             &
     &  ,Stashwork(pt401),Stashwork(pt402),Stashwork(pt403)             &
     &  ,Stashwork(pt404),Stashwork(pt405),Stashwork(pt406)             &
     &  ,Stashwork(pt407),Stashwork(pt408),Stashwork(pt409)             &
     &  ,Stashwork(pt410),Stashwork(pt411),Stashwork(pt412)             &
     &  ,Stashwork(pt413),Stashwork(pt414),Stashwork(pt415)             &
     &  ,Stashwork(pt416),Stashwork(pt417),Stashwork(pt418)             &
     &  ,Stashwork(pt419),stashwork(pt420),Stashwork(pt421)             &
     &  ,Stashwork(pt422),stashwork(pt423),Stashwork(pt424)             &
     &  ,Stashwork(pt425),stashwork(pt426),Stashwork(pt427)             &
     &  ,Stashwork(pt428),stashwork(pt429),Stashwork(pt430)             &
     &  ,Stashwork(pt431),stashwork(pt432),Stashwork(pt433)             &
     &  ,Stashwork(pt434),stashwork(pt435),Stashwork(pt436)             &
     &  ,Stashwork(pt437),stashwork(pt438),Stashwork(pt439)             &
     &  ,Stashwork(pt440),Stashwork(pt441)                              &

     &  ,Stashwork(pt451),Stashwork(pt452),Stashwork(pt453)             &
     &  ,Stashwork(pt454)  )


      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EOT_DIAG',4)
      END IF


      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STASH   ',3)
      END IF

      IF(ICODE /= 0) THEN
        RETURN
      ENDIF

! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,30,STASHWORK,                                &
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
     &           ICODE,CMESSAGE)

      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STASH   ',4)
      END IF

  999 CONTINUE
      RETURN
      END SUBROUTINE ST_DIAG3
