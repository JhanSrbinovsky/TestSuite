
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_geto2a(                                         &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
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
  &  GET_STEP,cmessage)


  USE mod_prism
  USE oasis3_atm_data_mod

  USE auscom_cpl_data_mod
  USE dump_received

  IMPLICIT NONE

  !
  ! Description:
  ! Receive incoming coupling data from the ocean (NEMO) to be used
  ! to drive the atmosphere. 
  !
  ! Author: R. Hill
  !
  ! Current Code Owner : R. Hill
  !
  !=======================================================================

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
! DECOMPTP comdeck
!
! Description
!
! Magic numbers indicating decomposition types.
! These numbers are used to index the arrays defined in the
! DECOMPDB comdeck, and are required as an argument to
! the CHANGE_DECOMPOSITION subroutine.
!
! Current code owner : P.Burton
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 4.2       19/08/96  Original code.   P.Burton
! 4.3       17/02/97  Added new ocean decomposition decomp_nowrap_ocean
!                     which does not contain extra wrap points at
!                     start and end of row.                  P.Burton
! 5.5       04/08/00  Modification for parallelisation of WAM
!                   Author:  Bob Carruthers, Cray UK Inc(D.Holmes-Bell)

! Magic Numbers indicating decomposition types

      INTEGER                                                           &
     &  max_decomps                                                     &
                               ! maximum number of decompositions
     &, decomp_unset                                                    &
                               ! no decomposition selected
     &, decomp_standard_atmos                                           &
                               ! standard 2D atmosphere
!                              ! decomposition
     &, decomp_standard_ocean                                           &
                               ! standard 1D ocean decomposition
     &, decomp_nowrap_ocean                                             &
                               ! 1D ocean without extra wrap-around
!                              ! points at ends of each row
     &, decomp_smexe                                                    &
     &, decomp_standard_wave   ! standard 1D WAM Wave Model
                               ! decomposition

      PARAMETER (                                                       &
     &  max_decomps=5                                                   &
     &, decomp_unset=-1                                                 &
     &, decomp_standard_atmos=1                                         &
     &, decomp_standard_ocean=2                                         &
     &, decomp_nowrap_ocean=3                                           &
     &, decomp_smexe=4                                                  &
     &, decomp_standard_wave=5)

! End of DECOMPTP comdeck
! DECOMPDB comdeck
!
! Description:
!
! DECOMPDB comdeck (Decomposition Database) contains information
! describing the various decompositions used by the mpp-UM
! The CHANGE_DECOMPOSITION subroutine can be used to select
! a particular decomposition (which copies the appropriate
! decomposition information into the PARVARS common block).
!
! Requires comdeck PARVARS to be *CALLed before it.
!
! Current code owner : P.Burton
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 4.2       19/08/96  Original code.   P.Burton
! 5.0       12/04/99  - added dimension Nfld_max to decomp_db_glsize
!                     - added dimension NHalo_max to decomp_db_halosize
!                     - added dimension NHalo_max to decomp_db_g_lasize
!                     - added dimension Nfld_max to
!                       decomp_db_g_lasize
!                     - replace blsizep/u variables by blsize with new
!                       Nfld_max dimension
!                                                            P.Burton
! 5.1       27/01/00  - changed g_pe_index to g_pe_index_EW
!                     - added g_pe_index_NS
!                                                  P.Burton
! 5.3       27/01/00  - Moved data statement into blkdata. A Van der Wal
! 5.3       14/09/01  Added model_domain variable.  P.Burton
! 5.5       30/01/03  Generalised datastart. P.Selwood.

! Common blocks containing information about each decomposition
! (For description of variables see the PARVARS comdeck)

      INTEGER :: decomp_db_bound(Ndim_max,max_decomps)                  &
     &, decomp_db_sb_model_domain(max_decomps)
      INTEGER :: decomp_db_glsize(Ndim_max,Nfld_max,max_decomps)
      INTEGER :: decomp_db_gridsize(Ndim_max,max_decomps)
      INTEGER :: decomp_db_g_lasize(Ndim_max,Nfld_max,                  &
     &                     NHalo_max,0:maxproc,max_decomps)
      INTEGER :: decomp_db_g_blsize(Ndim_max,Nfld_max,0:maxproc,        &
     &                          max_decomps)
      INTEGER :: decomp_db_g_datastart(Ndim_max,0:maxproc,max_decomps)
      INTEGER :: decomp_db_g_datastart_f(Ndim_max,Nfld_max,0:maxproc,   &
     &                                   max_decomps)
      INTEGER :: decomp_db_g_gridpos(Ndim_max,0:maxproc,max_decomps)
      INTEGER :: decomp_db_g_pe_index_EW(1-Max_Halo_Size:               &
     &  ROW_LENGTH_MAX+Max_Halo_Size, max_decomps)
      INTEGER :: decomp_db_g_pe_index_NS(1-Max_Halo_Size:               &
     &  ROWS_MAX+Max_Halo_Size, max_decomps)
      INTEGER :: decomp_db_halosize(Ndim_max,NHalo_max,max_decomps)
      INTEGER :: decomp_db_neighbour(4,max_decomps)
      INTEGER :: decomp_db_first_comp_pe(max_decomps)
      INTEGER :: decomp_db_last_comp_pe(max_decomps)
      INTEGER :: decomp_db_nproc(max_decomps)
      INTEGER :: decomp_db_gc_proc_row_group(max_decomps)
      INTEGER :: decomp_db_gc_proc_col_group(max_decomps)
      INTEGER :: decomp_db_gc_all_proc_group(max_decomps)

      ! indicates if a decomposition has been initialised

      LOGICAL :: decomp_db_set(max_decomps)

      COMMON /DECOMP_DATABASE/                                          &
     &  decomp_db_bound , decomp_db_sb_model_domain , decomp_db_glsize  &
     & ,decomp_db_g_lasize , decomp_db_gridsize,                        &
     &  decomp_db_g_blsize,                                             &
     &  decomp_db_g_datastart , decomp_db_g_datastart_f,                &
     &  decomp_db_g_gridpos,                                            &
     &  decomp_db_g_pe_index_EW,decomp_db_g_pe_index_NS,                &
     &  decomp_db_halosize , decomp_db_neighbour,                       &
     &  decomp_db_first_comp_pe , decomp_db_last_comp_pe,               &
     &  decomp_db_nproc,                                                &
     &  decomp_db_gc_proc_row_group , decomp_db_gc_proc_col_group,      &
     &  decomp_db_gc_all_proc_group,                                    &
     &  decomp_db_set

! End of DECOMPDB comdeck
!====================== COMDECK ATM_LSM ========================
! Description:
!   This comdeck contains a COMMON block which contains the
!   atmosphere land sea mask - both the full field, and the
!   local subdomain on this processor.
!   This data is required for various compression/decompression
!   algorithms.
!
!   Requires AMAXSIZE comdeck to be called first for Max2DFieldSize
!
! History:
!   Model    Date     Modification history
!   version
!   4.5      12/01/98 New comdeck created.                P.Burton
!   5.2      15/01/01 Only define common block in atmos model. R. Hill
!   5.4      06/09/02 SX now uses this deck - add def UTILIO.
!                                                    E.Leung
!   5.5      17/02/03 Upgrade Wave model from 4.1 to 5.5 D.Holmes-Bell
!   6.2      23/11/05 Removed all references to the wavemodel.
!                     T.Edwards
!

      LOGICAL                                                           &
!  Full-grid land-sea mask:
     &  atmos_landmask(Max2DFieldSize)                                  &
! Local subdomain area land-sea mask:
     &, atmos_landmask_local(Max2DFieldSize)

      INTEGER atmos_number_of_landpts ! total number of land points

      COMMON /Atmos_LSM_Common/                                         &
     &  atmos_landmask                                                  &
     &, atmos_landmask_local                                            &
     &, atmos_number_of_landpts


! End of comdeck ATM_LSM

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
! TYPPTRA
! include TYPSIZE first
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add pointers to new slab prognostics. J F Thomson.
!  3.5   19/05/95  Remove pointers JK1,JK2,JEXPK1,JEXPK2,JKDA,JKDF
!                  and JRHCRIT. D. Robinson
!  4.1   13/10/95  Add pointers for new soil moisture fraction
!                  prognostics,canopy conductance and
!                  vegetation J.Smith
!  4.1   26/04/96  Add pointers for Sulphur Cycle (12)   MJWoodage
!  4.3    18/3/97  Add pointers for HadCM2 sulphate loading patterns
!                                                   William Ingram
!  4.4   05/09/97  Add pointer for net energy flux prognostic
!                  S.D.Mullerworth
!  4.4   05/08/97  Add pointer for conv. cld amt on model levs. JMG
!  4.4  10/09/97   Added pointers for snow grain size and snow soot
!                  content used in prognostic snow albedo scheme
!                                                        R. Essery
!  4.4    16/9/97  Add pointers for new vegetation and land surface
!                  prognostics.                       Richard Betts
!  4.5    1/07/98  Add pointer for ocean CO2 flux. C.D.Jones
!  4.5    19/01/98 Replace JVEG_FLDS and JSOIL_FLDS with
!                  individual pointers. D. Robinson
!  4.5    04/03/98 Add 2 pointers for NH3 in S Cycle     M. Woodage
!                  Add 5 pointers for soot               M. Woodage
!                  Add pointer SO2_HILEM to CARGPT_ATMOS M. Woodage
!  4.5    08/05/98 Add 16 pointers for User Anc.         D. Goddard
!  4.5    15/07/98 Add pointers for new 3D CO2 array.    C.D.Jones
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  5.1    25/02/00 Add new LBC pointers and labelled old JRIM
!                  pointers ready for removal.            P.Burton
!  5.1    28/04/00 Add extra level to p/exner on rho levels
!                                                 Andy Malcolm
!  5.1    29/02/00 Added JNET_MFLUX pointer for net moisture flux
!  5.1    03/05/00 Add JTSTAR_ANOM pointer to CARGPT_ATMOS. D Robinson
!  5.2    18/10/00 Inserted lines for slab model fields   K.D.Williams
!  5.1    31/08/00 Added JOROG_LBC pointer
!                  Move JEXNER_LBC to top of list
!                  Remove levels from the LBC J_pointers    P.Burton
!  5.2    25/09/00 Clear out diagnostic RHcrit legacy code. A.C.Bushell
!
!  5.2    13/09/00 Remove JSOOT, JSO4, JH2SO4 pointers. P.Selwood.
!  5.2    15/11/00 Re-Introduce pointers for MOSES 2 and Triffid
!                  M. Best.
!  5.2    09/03/01 Introduce extra pointers in level dependent constants
!                  for height definitions. R Rawlins
!  5.3    04/10/01 Added new slab prognostics             K.Williams
!  5.3       06/01 Introduce extra pointers needed for coastal
!                  tiling.                              Nic Gedney
!  5.3    19/06/01 Added JTPPSOZONE (tropopauase-based ozone) Dave Tan
!  5.3    01/11/01 Add J_CO2_EMITS and J_CO2FLUX to common block
!  5.4    02/09/02 Added new pointer JSNODEP_SEA        K.Williams
!  5.4    28/08/02 Add JCATCH_SNOW and JSNOW_GRND to common block
!  5.5    05/02/03 Add pointers for biomass aerosol scheme. P Davison
!  5.5    05/11/02 Large-scale hydrology scheme: Add pointers for
!                  gamma function, water table depth,
!                  surface saturated fraction and wetland fraction.
!                                                        N. Gedney
!  5.5    24/02/03 Add JICE_FRACT_CAT, JICE_THICK_CAT, JTI_CAT &
!                  JSNODEP_SEA_CAT to common block      J.Ridley
!  5.5    03/02/03 Add JQCF2,JQRAIN,JQGRAUP and LBC variables. R.Forbes
!  5.5    26/02/03 River routing support. P.Selwood.
!  6.1    07/09/04 Add pointers for STOCHEM fields.     C. Johnson
!  6.0    30/07/03 Add cloud fraction lbc variables. Damian Wilson
!  6.0    12/08/03 removes JRIV_GRIDBOX. C.Bunton.
!  6.0    12/09/03 Extra pointers for river routing 2A. V.A.Bell
!  6.1    07/04/04 Add pointer for DMS concentration.   A. Jones
!  6.2    24/02/06 Add pointer for DMS flux from ocean model    J.Gunson
!   6.2   21/2/06  Declare new pointer to
!                  re-route outflow from inland basins to soil moisture
!                  P. Falloon
!  6.2    01/03/06 Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2    19/01/06 Add variables required for surface radiative
!                  forcing. Needed under versions 3C and 3Z of
!                  the radiation code.            (J.-C.Thelen)
!  6.2    01/10/05 Add JMURK_LBC variables.   R.M.Forbes
!  6.2    14/11/05 Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2    01/03/06 Add pointer for RothC pools and fluxes.   C.D. Jones
! Pointers for ATMOSPHERE model variables. Configuration dependent.
! Addresses in D1 array of model variables held in primary and
! secondary space.

      ! 1: Array Variables (dimensions are resolution dependent.)

      ! 1.1: Data variables stored in primary space.
      INTEGER :: JU(MODEL_LEVELS)      ! u component of wind
      INTEGER :: JV(MODEL_LEVELS)      ! v component of wind
      INTEGER :: JW(0:MODEL_LEVELS)    ! w component of wind
      INTEGER :: JRHO(MODEL_LEVELS)    ! Density
      INTEGER :: JTHETA(MODEL_LEVELS)  ! Potential temperature
      INTEGER :: JQ(WET_LEVELS)        ! Specific humidity
      INTEGER :: JQCL(WET_LEVELS)      ! qcl
      INTEGER :: JQCF(WET_LEVELS)      ! qcf
      INTEGER :: JQCF2(WET_LEVELS)     ! second ice
      INTEGER :: JQRAIN(WET_LEVELS)    ! rain
      INTEGER :: JQGRAUP(WET_LEVELS)   ! graupel

      ! Exner pressure on rho levels
      INTEGER :: JEXNER_RHO_LEVELS(MODEL_LEVELS+1)

      INTEGER :: JU_ADV(MODEL_LEVELS)   ! Advective u component of wind
      INTEGER :: JV_ADV(MODEL_LEVELS)   ! Advective v component of wind
      INTEGER :: JW_ADV(0:MODEL_LEVELS) ! Advective w component of wind

      ! 1.2: Data variables stored in secondary space.
      INTEGER :: JP(MODEL_LEVELS+1)                  ! Pressure on rho le

      ! Pressure on theta levels
      INTEGER :: JP_THETA_LEVELS(MODEL_LEVELS)

      ! Exner pressure on theta levels
      INTEGER :: JEXNER_THETA_LEVELS(MODEL_LEVELS)

      ! 1.3: Cloud Fields
      INTEGER :: JCCW_RAD(WET_LEVELS)            ! CCW profile to radiation
      INTEGER :: JCCA(N_CCA_LEV)                 ! Convective cloud amount
      INTEGER :: JCF_AREA(WET_LEVELS)            ! Area Cloud Fraction
      INTEGER :: JCF_BULK(WET_LEVELS)            ! Bulk Cloud Fraction
      INTEGER :: JCF_LIQUID(WET_LEVELS)          ! Liquid cloud fraction
      INTEGER :: JCF_FROZEN(WET_LEVELS)          ! Frozen cloud fraction

      ! 1.4: Soil Ancillary fields
      INTEGER :: J_DEEP_SOIL_TEMP(ST_LEVELS)   ! Deep soil temperature

      INTEGER :: JSMCL(SM_LEVELS)   ! Soil moisture content in layers
      INTEGER :: JSTHU(SM_LEVELS)   ! Unfrozen soil moisture fraction
      INTEGER :: JSTHF(SM_LEVELS)   ! Frozen soil moisture fraction

      ! 1.5: Radiation Increments
      INTEGER :: JSW_INCS(0:MODEL_LEVELS+1)    ! SW radiation increments
      INTEGER :: JLW_INCS(0:MODEL_LEVELS)      ! LW radiation increments
! PAR radiation increment
      INTEGER :: JDIRPAR

      ! 1.6: Ozone
      INTEGER :: JOZONE(OZONE_LEVELS)          ! Ozone
!  tropopause-based ozone
      INTEGER :: JTPPSOZONE(tpps_ozone_levels)
      ! 1.7: Tracer and aerosol fields
      INTEGER :: JTRACER(TR_LEVELS,TR_VARS+1)  ! Tracers
      INTEGER :: JTR_UKCA(TR_LEVELS,TR_UKCA+1)  ! UKCA Tracers
      INTEGER :: JMURK_SOURCE(MODEL_LEVELS)    ! multi-level murk source
      INTEGER :: JMURK(MODEL_LEVELS)           ! multi-level murk concent
      INTEGER :: JDUST_DIV1(MODEL_LEVELS)      ! dust mmr, division 1
      INTEGER :: JDUST_DIV2(MODEL_LEVELS)      ! dust mmr, division 2
      INTEGER :: JDUST_DIV3(MODEL_LEVELS)      ! dust mmr, division 3
      INTEGER :: JDUST_DIV4(MODEL_LEVELS)      ! dust mmr, division 4
      INTEGER :: JDUST_DIV5(MODEL_LEVELS)      ! dust mmr, division 5
      INTEGER :: JDUST_DIV6(MODEL_LEVELS)      ! dust mmr, division 6
      INTEGER :: JSO2(MODEL_LEVELS)            ! sulphur dioxide gas
      INTEGER :: JDMS(MODEL_LEVELS)            ! dimethyl sulphide gas
      INTEGER :: JSO4_AITKEN(MODEL_LEVELS)     ! Aitken mode sulphate aer
      INTEGER :: JSO4_ACCU(MODEL_LEVELS)       ! accumulation mode sulpha
      INTEGER :: JSO4_DISS(MODEL_LEVELS)       ! dissloved  sulphate aero
      INTEGER :: JH2O2(MODEL_LEVELS)           ! hydrogen peroxide mmr
      INTEGER :: JNH3(MODEL_LEVELS)            ! ammonia gas mmr
      INTEGER :: JSOOT_NEW(MODEL_LEVELS)       ! fresh soot mmr
      INTEGER :: JSOOT_AGD(MODEL_LEVELS)       ! aged soot mmr
      INTEGER :: JSOOT_CLD(MODEL_LEVELS)       ! soot in cloud mmr
      INTEGER :: JBMASS_NEW(MODEL_LEVELS)      ! fresh biomass mmr
      INTEGER :: JBMASS_AGD(MODEL_LEVELS)      ! aged biomass mmr
      INTEGER :: JBMASS_CLD(MODEL_LEVELS)      ! cloud biomass mmr
      INTEGER :: JOCFF_NEW(MODEL_LEVELS)       ! fresh OCFF mmr
      INTEGER :: JOCFF_AGD(MODEL_LEVELS)       ! aged OCFF mmr
      INTEGER :: JOCFF_CLD(MODEL_LEVELS)       ! OCFF in cloud mmr
      INTEGER :: JSO2_NATEM(MODEL_LEVELS)      ! natural SO2 emissions
      INTEGER :: JOH(MODEL_LEVELS)             ! hydroxyl radical ancilla
      INTEGER :: JHO2(MODEL_LEVELS)            ! hydrogen dioxide ancilla
      INTEGER :: JH2O2_LIMIT(MODEL_LEVELS)     ! limiting H2O2 ancillary
      INTEGER :: JO3_CHEM(MODEL_LEVELS)        ! ozone for chemistry anci
      INTEGER :: JHadCM2_SO4(2)                ! HadCM2 sulphate loading
      ! JHadCM2_SO4: Should really be NSULPHAT (= 2) but to use NSULPAT,
      !              must include CSENARIO before every include of TYPPTR
      INTEGER :: JCO2(MODEL_LEVELS)            ! 3D CO2 FIELD
      INTEGER :: JCH4_STOCH(MODEL_LEVELS)      ! STOCHEM CH4
      INTEGER :: JO3_STOCH(MODEL_LEVELS)       ! STOCHEM O3
      INTEGER :: JOH_UKCA(MODEL_LEVELS)        ! OH MMR from UKCA
      INTEGER :: JHO2_UKCA(MODEL_LEVELS)       ! HO2 MMR from UKCA
      INTEGER :: JH2O2_UKCA(MODEL_LEVELS)      ! H2O2 MMR from UKCA
      INTEGER :: JO3_UKCA(MODEL_LEVELS)        ! O3 MMR from UKCA
      INTEGER :: JOZONE_TRACER(MODEL_LEVELS)   ! Prognostic O3 Tracer(Cariol)
      INTEGER :: JO3_PROD_LOSS(MODEL_LEVELS)   ! Cariol O3 Prod-Loss (P-L)
      INTEGER :: JO3_P_L_VMR(MODEL_LEVELS)     ! Cariol O3 P-L wrt VMR 
      INTEGER :: JO3_VMR(MODEL_LEVELS)         ! Cariol O3 Vol Mix Ratio-VMR
      INTEGER :: JO3_P_L_TEMP(MODEL_LEVELS)    ! Cariol O3 P-L wrt temp 
      INTEGER :: JO3_TEMP(MODEL_LEVELS)        ! Cariol O3 temp  
      INTEGER :: JO3_P_L_COLO3(MODEL_LEVELS)   ! Cariol O3 P-L wrt colO3 
      INTEGER :: JO3_COLO3(MODEL_LEVELS)       ! Cariol O3 column (colO3) 
      INTEGER :: JARCLBIOG_BG(MODEL_LEVELS)    ! Biogenic aerosol climatology 
      INTEGER :: JARCLBIOM_FR(MODEL_LEVELS)    ! Biomass burning (fresh) aerosol clim
      INTEGER :: JARCLBIOM_AG(MODEL_LEVELS)    ! Biomass burning (aged) aerosol clim
      INTEGER :: JARCLBIOM_IC(MODEL_LEVELS)    ! Biomass burning (in-cloud) aerosol clim
      INTEGER :: JARCLBLCK_FR(MODEL_LEVELS)    ! Black carbon (fresh) aerosol clim
      INTEGER :: JARCLBLCK_AG(MODEL_LEVELS)    ! Black carbon (aged) aerosol clim
      INTEGER :: JARCLSSLT_FI(MODEL_LEVELS)    ! Sea salt (film mode) aerosol clim 
      INTEGER :: JARCLSSLT_JT(MODEL_LEVELS)    ! Sea salt (jet mode) aerosol clim
      INTEGER :: JARCLSULP_AC(MODEL_LEVELS)    ! Sulphate (accumulation mode) aero clim
      INTEGER :: JARCLSULP_AK(MODEL_LEVELS)    ! Sulphate (Aitken mode) aerosol clim 
      INTEGER :: JARCLSULP_DI(MODEL_LEVELS)    ! Sulphate (dissolved) aerosol clim
      INTEGER :: JARCLDUST_B1(MODEL_LEVELS)    ! Dust (bin 1) aerosol climatology 
      INTEGER :: JARCLDUST_B2(MODEL_LEVELS)    ! Dust (bin 2) aerosol climatology 
      INTEGER :: JARCLDUST_B3(MODEL_LEVELS)    ! Dust (bin 3) aerosol climatology 
      INTEGER :: JARCLDUST_B4(MODEL_LEVELS)    ! Dust (bin 4) aerosol climatology 
      INTEGER :: JARCLDUST_B5(MODEL_LEVELS)    ! Dust (bin 5) aerosol climatology 
      INTEGER :: JARCLDUST_B6(MODEL_LEVELS)    ! Dust (bin 6) aerosol climatology 
      INTEGER :: JARCLOCFF_FR(MODEL_LEVELS)    ! Org carbon fossil fuel (fresh) aero clim
      INTEGER :: JARCLOCFF_AG(MODEL_LEVELS)    ! Org carbon fossil fuel (aged) aero clim
      INTEGER :: JARCLOCFF_IC(MODEL_LEVELS)    ! Org carbon fossil fuel (in-cloud) aero clim
      INTEGER :: JARCLDLTA_DL(MODEL_LEVELS)    ! Delta aerosol climatology
!

! 1.8: Multi-level user ancillary fields
      INTEGER :: JUSER_MULT1(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT2(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT3(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT4(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT5(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT6(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT7(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT8(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT9(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT10(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT11(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT12(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT13(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT14(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT15(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT16(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT17(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT18(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT19(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT20(MODEL_LEVELS)    ! multi-level user ancilla

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
      INTEGER :: JOROG_LBC                       ! Orography LBC
      INTEGER :: JU_LBC                          ! U LBC
      INTEGER :: JV_LBC                          ! V LBC
      INTEGER :: JW_LBC                          ! W LBC
      INTEGER :: JRHO_LBC                        ! RHO LBC
      INTEGER :: JTHETA_LBC                      ! Theta LBC
      INTEGER :: JQ_LBC                          ! Q LBC
      INTEGER :: JQCL_LBC                        ! QCL LBC
      INTEGER :: JQCF_LBC                        ! QCF LBC
      INTEGER :: JQCF2_LBC                       ! 2nd Ice LBC
      INTEGER :: JQRAIN_LBC                      ! Rain LBC
      INTEGER :: JQGRAUP_LBC                     ! Graupel LBC
      INTEGER :: JCF_BULK_LBC                    ! CF_BULK LBC
      INTEGER :: JCF_LIQUID_LBC                  ! CF_LIQUID_LBC
      INTEGER :: JCF_FROZEN_LBC                  ! CF_FROZEN_LBC
      INTEGER :: JEXNER_LBC                      ! Exner LBC
      INTEGER :: JU_ADV_LBC                      ! U_ADV LBC
      INTEGER :: JV_ADV_LBC                      ! V_ADV LBC
      INTEGER :: JW_ADV_LBC                      ! W_ADV LBC
      INTEGER :: JMURK_LBC                       ! Murk aerosol LBC
      INTEGER :: JTRACER_LBC(TR_VARS+1)          ! Tracer LBCs
      ! Lateral Boundary Condition tendencies
      INTEGER :: JU_LBC_TEND                     ! U LBC  tendencies
      INTEGER :: JV_LBC_TEND                     ! V LBC tendencies
      INTEGER :: JW_LBC_TEND                     ! W LBC tendencies
      INTEGER :: JRHO_LBC_TEND                   ! RHO LBC tendencies
      INTEGER :: JTHETA_LBC_TEND                 ! Theta LBC tendencies
      INTEGER :: JQ_LBC_TEND                     ! Q LBC tendencies
      INTEGER :: JQCL_LBC_TEND                   ! QCL LBC tendencies
      INTEGER :: JQCF_LBC_TEND                   ! QCF LBC tendencies
      INTEGER :: JQCF2_LBC_TEND                  ! 2nd Ice
      INTEGER :: JQRAIN_LBC_TEND                 ! Rain LBC tendencies
      INTEGER :: JQGRAUP_LBC_TEND                ! Graupel
      INTEGER :: JCF_BULK_LBC_TEND               ! CF_BULK LBC tend'cies
      INTEGER :: JCF_LIQUID_LBC_TEND             ! CF_LIQUID_LBC t'cies
      INTEGER :: JCF_FROZEN_LBC_TEND             ! CF_FROZEN_LBC t'cies
      INTEGER :: JEXNER_LBC_TEND                 ! Exner LBC tendencies
      INTEGER :: JU_ADV_LBC_TEND                 ! U_ADV LBC tendencies
      INTEGER :: JV_ADV_LBC_TEND                 ! V_ADV LBC tendencies
      INTEGER :: JW_ADV_LBC_TEND                 ! W_ADV LBCtendencies
      INTEGER :: JMURK_LBC_TEND                  ! Murk aerosol LBC tend
      INTEGER :: JTRACER_LBC_TEND(TR_VARS+1)     ! Tracer LBC tendencies

      ! 2: Scalar Variables

      ! 2.1: Data variables stored in primary space.
      INTEGER :: JTSTAR         ! Surface temperature
      INTEGER :: JLAND          ! Land sea mask
      INTEGER :: JTSTAR_ANOM    ! Surface temperature anolomy
!   2.15: Fields for coastal tiling
      INTEGER :: JFRAC_LAND  ! Land fraction in grid box
      INTEGER :: JTSTAR_LAND ! Land surface temperature
      INTEGER :: JTSTAR_SEA  ! Sea surface temperature
      INTEGER :: JTSTAR_SICE ! Sea-ice surface temperature
! Set pointers for sea-ice and land albedos
      INTEGER :: JSICE_ALB                   ! Sea-ice albedo
      INTEGER :: JLAND_ALB                   ! Mean land albedo

      ! 2.2: Data variables stored in secondary space.

      INTEGER :: JPSTAR          ! Surface pressure

      ! 2.3: Cloud fields
      INTEGER :: JLCBASE         ! Lowest Convective cloud base
      INTEGER :: JCCB            ! Convective cloud base
      INTEGER :: JCCT            ! Convective cloud top

      INTEGER :: JCCLWP          ! Convective cloud liquid water path

      ! 2.4: Boundary layer fields

      INTEGER :: JZH                         ! Boundary layer depth

      ! Standard deviation of turbulent fluctuations of layer 1
                                    ! temperature
      INTEGER :: JT1_SD

      ! Standard deviation of turbulent fluctuations of layer 1 humidity
      INTEGER :: JQ1_SD

      ! Number of model levels in the  turbulently mixed layer
      INTEGER :: JNTML

      ! Top level for turb mixing in any decoupled Sc layer
      INTEGER :: JNTDSC

      ! Bottom level for turb mixing in any decoupled Sc layer
      INTEGER :: JNBDSC

      INTEGER :: JCUMULUS      ! Boundary layer convection flag

      ! 2.4: Soil Ancillary fields

      INTEGER :: JSAT_SOILW_SUCTION          ! Saturated soil water sucti
      INTEGER :: JTHERM_CAP     ! Thermal capacity
      INTEGER :: JTHERM_COND    ! Thermal conductivity
      INTEGER :: JVOL_SMC_CRIT  ! Vol smc at critical point
      INTEGER :: JVOL_SMC_WILT  ! Vol smc at wilting point
      INTEGER :: JVOL_SMC_SAT   ! Vol smc at saturation
      INTEGER :: JSAT_SOIL_COND ! Saturated soil conductivity
      INTEGER :: JCLAPP_HORN    ! Clapp-Hornberger B coefficient

      ! 2.5: Vegetation Ancillary fields

      INTEGER :: JCANOPY_WATER  ! Canopy Water
      INTEGER :: JSURF_CAP      ! Surface Capacity
      INTEGER :: JSURF_RESIST   ! Surface Resistance
      INTEGER :: JROOT_DEPTH    ! Root depth
      INTEGER :: JINFILT        ! Infiltration factor
      INTEGER :: JVEG_FRAC      ! Vegetation fraction
      INTEGER :: JLAI           ! Leaf area index
      INTEGER :: JCANHT         ! Canopy height
      INTEGER :: JZ0            ! Vegetative Roughness length
      INTEGER :: JSFA           ! Snow free albedo
      INTEGER :: JMDSA          ! Cold deep snow albedo
      INTEGER :: JGS            ! Gridbox mean canopy conductance

      ! 2.6: Orographic Ancillary fields

      INTEGER :: JOROG          ! Orographic height
      INTEGER :: JOROG_SD       ! Standard Deviation of orography
      INTEGER :: JOROG_SIL      ! Silhouette area of orography
      INTEGER :: JOROG_HO2      ! Peak to trough height/(2*sqrt2)
      INTEGER :: JOROG_GRAD_X   ! Orographic gradient x
      INTEGER :: JOROG_GRAD_Y   ! Orographic gradient y
      INTEGER :: JOROG_GRAD_XX  ! Orographic gradient xx
      INTEGER :: JOROG_GRAD_XY  ! Orographic gradient xy
      INTEGER :: JOROG_GRAD_YY  ! Orographic gradient yy

      ! 2.7: Sea/Sea Ice

      INTEGER :: JU_SEA         ! Surface current (u component)
      INTEGER :: JV_SEA         ! Surface current (v component)
      INTEGER :: JU_0_P         ! Surace  current (u) on p-grid
      INTEGER :: JV_0_P         ! Surface current (v) on p-grid
      INTEGER :: JICE_FRACTION  ! Sea ice fraction
      INTEGER :: JICE_THICKNESS ! Sea ice depth
      INTEGER :: JTI            ! Sea ice temperature
      INTEGER :: JICE_FRACT_CAT ! Sea ice fraction on catagories
      INTEGER :: JICE_THICK_CAT ! Sea ice thickness on catagories
      INTEGER :: JTI_CAT        ! Sea ice temperature on catagories

      ! 2.8: Snow

      INTEGER :: JSNODEP        ! Snow depth on land
      INTEGER :: JSNODEP_SEA    ! Snow depth on sea ice
      INTEGER :: JSNODEP_SEA_CAT ! Snow depth on sea ice catagories
      INTEGER :: JCATCH_SNOW    ! Coniferous canopy
!                               ! snow capacity
      INTEGER :: JSNOW_GRND     ! Snow below canopy
      INTEGER :: JSNSOOT        ! Snow soot content

! 2.9: aerosol emission fields, including mineral dust parent soil props

      INTEGER :: JSOIL_CLAY                    ! soil clay fraction
      INTEGER :: JSOIL_SILT                    ! soil silt fraction
      INTEGER :: JSOIL_SAND                    ! soil sand fraction
      INTEGER :: JDUST_MREL1                   ! soil rel mass, div 1
      INTEGER :: JDUST_MREL2                   ! soil rel mass, div 2
      INTEGER :: JDUST_MREL3                   ! soil rel mass, div 3
      INTEGER :: JDUST_MREL4                   ! soil rel mass, div 4
      INTEGER :: JDUST_MREL5                   ! soil rel mass, div 5
      INTEGER :: JDUST_MREL6                   ! soil rel mass, div 6


      INTEGER :: JSO2_EM        ! sulphur dioxide emission
      INTEGER :: JDMS_EM        ! dimethyl sulphide emission
      INTEGER :: JSO2_HILEM     ! high level SO2 emissions
      INTEGER :: JNH3_EM        ! ammonia gas surface emiss
      INTEGER :: JSOOT_EM       ! fresh soot surface emissions
      INTEGER :: JSOOT_HILEM    ! fresh soot high lev emissions
      INTEGER :: JBMASS_EM       ! fresh bmass surface emissions
      INTEGER :: JBMASS_HILEM    ! fresh bmass high lev emissions
      INTEGER :: JOCFF_EM        ! fresh OCFF surface emissions
      INTEGER :: JOCFF_HILEM     ! fresh OCFF high-level emissions
      INTEGER :: JDMS_CONC       ! seawater dimethyl sulphide conc.
      INTEGER :: JDMS_OFLUX      ! DMS flux from ocean model

      ! 2.10: User ancillary fields
      INTEGER :: JUSER_ANC1         ! user ancillary field 1
      INTEGER :: JUSER_ANC2                  ! user ancillary field 2
      INTEGER :: JUSER_ANC3                  ! user ancillary field 3
      INTEGER :: JUSER_ANC4                  ! user ancillary field 4
      INTEGER :: JUSER_ANC5                  ! user ancillary field 5
      INTEGER :: JUSER_ANC6                  ! user ancillary field 6
      INTEGER :: JUSER_ANC7                  ! user ancillary field 7
      INTEGER :: JUSER_ANC8                  ! user ancillary field 8
      INTEGER :: JUSER_ANC9                  ! user ancillary field 9
      INTEGER :: JUSER_ANC10                 !user ancillary field 10
      INTEGER :: JUSER_ANC11                 ! user ancillary field 11
      INTEGER :: JUSER_ANC12                 ! user ancillary field 12
      INTEGER :: JUSER_ANC13                 ! user ancillary field 13
      INTEGER :: JUSER_ANC14                 ! user ancillary field 14
      INTEGER :: JUSER_ANC15                 ! user ancillary field 15
      INTEGER :: JUSER_ANC16                 ! user ancillary field 16
      INTEGER :: JUSER_ANC17                 ! user ancillary field 17
      INTEGER :: JUSER_ANC18                 ! user ancillary field 18
      INTEGER :: JUSER_ANC19                 ! user ancillary field 19
      INTEGER :: JUSER_ANC20                 ! user ancillary field 20

      ! Tracer Fluxes - kdcorbin, 05/10
      INTEGER :: JTRACER_FLUX1
      INTEGER :: JTRACER_FLUX2
      INTEGER :: JTRACER_FLUX3
      INTEGER :: JTRACER_FLUX4
      INTEGER :: JTRACER_FLUX5
      INTEGER :: JTRACER_FLUX6
      INTEGER :: JTRACER_FLUX7
      INTEGER :: JTRACER_FLUX8
      INTEGER :: JTRACER_FLUX9
      INTEGER :: JTRACER_FLUX10
      INTEGER :: JTRACER_FLUX11
      INTEGER :: JTRACER_FLUX12
      INTEGER :: JTRACER_FLUX13
      INTEGER :: JTRACER_FLUX14
      INTEGER :: JTRACER_FLUX15
      INTEGER :: JTRACER_FLUX16
      INTEGER :: JTRACER_FLUX17
      INTEGER :: JTRACER_FLUX18
      INTEGER :: JTRACER_FLUX19
      INTEGER :: JTRACER_FLUX20  

      !   2.11: Store arrays for energy correction calculation
      INTEGER :: JNET_FLUX                   ! Net energy flux
      INTEGER :: JNET_MFLUX                  ! Net moisture flux

      !   2.12: Tiled Vegetation and Triffid fields
      INTEGER :: JFRAC_TYP        ! Fractions of surface type
      INTEGER :: JFRAC_CON1       ! Fractions of surface type
      INTEGER :: JFRAC_CON2       ! Fractions of surface type
      INTEGER :: JFRAC_CON3       ! Fractions of surface type
      INTEGER :: JFRAC_CON4       ! Fractions of surface type
      INTEGER :: JFRAC_CON5       ! Fractions of surface type
      INTEGER :: JFRAC_CON6       ! Fractions of surface type
      INTEGER :: JFRAC_CON7       ! Fractions of surface type
      INTEGER :: JFRAC_CON8       ! Fractions of surface type
      INTEGER :: JFRAC_CON9       ! Fractions of surface type
      INTEGER :: JLAI_PFT         ! LAI of plant functional types
      INTEGER :: JCANHT_PFT       ! Canopy hght of plant func types
      INTEGER :: JDISTURB         ! Disturbed fraction of vegetation
      INTEGER :: JSOIL_ALB        ! Snow-free albedo of bare soil
      INTEGER :: JSOIL_CARB       ! Soil carbon content
      INTEGER :: JSOIL_CARB1      ! Soil carbon content DPM
      INTEGER :: JSOIL_CARB2      ! Soil carbon content RPM
      INTEGER :: JSOIL_CARB3      ! Soil carbon content BIO
      INTEGER :: JSOIL_CARB4      ! Soil carbon content HUM
      INTEGER :: JNPP_PFT_ACC     ! Accumulated NPP on PFTs
      INTEGER :: JG_LF_PFT_ACC    ! Accum. leaf turnover rate PFTs
      INTEGER :: JG_PHLF_PFT_ACC  ! Accumulated phenological leaf
                                    ! turnover rate on PFTs
      INTEGER :: JRSP_W_PFT_ACC   ! Accum. wood respiration on PFTs
      INTEGER :: JRSP_S_ACC       ! Accumulated soil respiration
      INTEGER :: JRSP_S_ACC1      ! Accumulated soil respiration DPM
      INTEGER :: JRSP_S_ACC2      ! Accumulated soil respiration RPM
      INTEGER :: JRSP_S_ACC3      ! Accumulated soil respiration BIO
      INTEGER :: JRSP_S_ACC4      ! Accumulated soil respiration HUM
      INTEGER :: JCAN_WATER_TILE  ! Canopy water content on tiles
      INTEGER :: JCATCH_TILE      ! Canopy capacity on tiles
      INTEGER :: JINFIL_TILE      ! Max infiltration rate on tiles
      INTEGER :: JRGRAIN_TILE     ! Snow grain size on tiles
      INTEGER :: JSNODEP_TILE     ! Snow depth on tiles
      INTEGER :: JTSTAR_TILE      ! Surface temperature on tiles
      INTEGER :: JZ0_TILE         ! Surface roughness on tiles
      INTEGER :: JDOLR            ! TOA - surface upward LW at
                                  ! radiation timestep
      INTEGER :: JLW_DOWN         ! Surface downward LW at
                                  ! radiation timestep
      INTEGER :: JSW_TILE         ! Surface net SW on land tiles at
                                  ! radiation timestep
      !   2.13: Slab Model
      INTEGER :: JTSLAB           ! Temperature of slab ocean.
      INTEGER :: JTCLIM           ! Climatological SST's
      INTEGER :: JHCLIM           ! Climatological ice depth
      INTEGER :: JCHEAT     ! Caryheat (from ice model to ocn)
      INTEGER :: JOIFLX     ! Ocean-Ice heat flux
      INTEGER :: JUICE      ! Zonal comp of ice velocity
      INTEGER :: JVICE      ! Meridional comp of ice velocity
      INTEGER :: JSIG11NE   !  Internal stresses for
      INTEGER :: JSIG11SE   !  EVP ice dynamics
      INTEGER :: JSIG11SW
      INTEGER :: JSIG11NW
      INTEGER :: JSIG12NE
      INTEGER :: JSIG12SE
      INTEGER :: JSIG12SW
      INTEGER :: JSIG12NW
      INTEGER :: JSIG22NE
      INTEGER :: JSIG22SE
      INTEGER :: JSIG22SW
      INTEGER :: JSIG22NW

!   2.14: Carbon cycle fields
      INTEGER :: J_CO2FLUX   ! Ocean CO2 flux (Kg CO2/m2/s1)
      INTEGER :: J_CO2_EMITS ! Surface CO2 emissions (Kg CO2/m2/s1)

!   2.15: Fields carried forward from previous version.
!         May not be required
      INTEGER :: JSURF_RESIST_NIT  ! Surface resistance on
                                    ! non-ice tiles
      INTEGER :: JROOT_DEPTH_NIT   ! Root depth on non-ice tiles
      INTEGER :: JZ0V_TYP          ! Vegetative Roughness length on
                                    ! tiles
      INTEGER :: JTSNOW            ! Snow surface layer temperature
      INTEGER :: JICE_EDGE
      INTEGER :: JOROG_TENDENCY    ! Orographic tendencies
      INTEGER :: JOROG_SD_TENDENCY ! Orographic variable tendency

      ! Pointers for ATMOSPHERE model constants. Scalars only.
      ! Addresses in level dependent constants array.
      INTEGER :: JETATHETA
      INTEGER :: JETARHO
      INTEGER :: JRHCRIT
      INTEGER :: JSOIL_THICKNESS
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      INTEGER :: Jzseak_theta ! zsea(k) on theta levels
      INTEGER :: JCk_theta    ! C(k)    on theta levels
      INTEGER :: Jzseak_rho   ! zsea(k) on rho levels
      INTEGER :: JCk_rho      ! C(k)    on rho levels
      ! Addresses in Row and Col  dependent constants array.
      INTEGER :: JLAMBDA_INPUT_P
      INTEGER :: JLAMBDA_INPUT_U
      INTEGER :: JPHI_INPUT_P
      INTEGER :: JPHI_INPUT_V
!   2.16: Fields for large-scale hydrology scheme.

      INTEGER :: JTI_MEAN          !Mean topographic index
      INTEGER :: JTI_SIG           !Standard dev. in topographic index
      INTEGER :: JFEXP             !Exponential decay in soil
!                                  ! saturated conductivity
      INTEGER :: JGAMMA_INT        !Integrated gamma function
      INTEGER :: JWATER_TABLE      !Water table depth
      INTEGER :: JFSFC_SAT         !Surface saturation fraction
      INTEGER :: JF_WETLAND        !Wetland fraction
      INTEGER :: JSTHZW           !Soil moist fract. in deep-zw layer.
      INTEGER :: JA_FSAT          !Fitting parameter for Fsat in LSH.
      INTEGER :: JC_FSAT          !Fitting parameter for Fsat in LSH.
      INTEGER :: JA_FWET          !Fitting parameter for Fwet in LSH.
      INTEGER :: JC_FWET          !Fitting parameter for Fwet in LSH.


!   2.17: Fields for River routing.
      INTEGER :: JRIV_SEQUENCE   ! River sequence
      INTEGER :: JRIV_DIRECTION  ! River direction
      INTEGER :: JRIV_STORAGE    ! River water storage
      INTEGER :: JTOT_SURFROFF   ! Accumulated surface runoff
      INTEGER :: JTOT_SUBROFF    !     "       sub-surface runoff
      INTEGER :: JRIV_INLANDATM       ! inland basin outflow
! Fields for grid-to-grid river routing (river routing 2A)
      INTEGER :: JRIV_IAREA      ! Drainage area
      INTEGER :: JRIV_SLOPE      ! Grid-cell slope
      INTEGER :: JRIV_FLOWOBS1   ! Initial values of flow
      INTEGER :: JRIV_INEXT      ! Flow direction (x)
      INTEGER :: JRIV_JNEXT      ! Flow direction (y)
      INTEGER :: JRIV_LAND       ! Land-type (land/river/sea)
      INTEGER :: JRIV_SUBSTORE   ! Subsurface storage
      INTEGER :: JRIV_SURFSTORE  ! Surface storage
      INTEGER :: JRIV_FLOWIN     ! Surface inflow
      INTEGER :: JRIV_BFLOWIN    ! Subsurface inflow

      INTEGER :: JC_SOLAR 
      INTEGER :: JC_BLUE 
      INTEGER :: JC_DOWN 
      INTEGER :: JC_LONGWAVE 
      INTEGER :: JC_TAUX 
      INTEGER :: JC_TAUY 
      INTEGER :: JC_WINDMIX 
      INTEGER :: JC_SENSIBLE 
      INTEGER :: JC_SUBLIM 
      INTEGER :: JC_EVAP 
      INTEGER :: JC_BOTMELTN 
      INTEGER :: JC_TOPMELTN 
      INTEGER :: JC_LSRAIN 
      INTEGER :: JC_LSSNOW 
      INTEGER :: JC_CVRAIN 
      INTEGER :: JC_CVSNOW 
      INTEGER :: JC_RIVEROUT 

      INTEGER :: JC_PRESS

!    Fields for CABLE
      INTEGER :: JTSOIL_TILE(SM_LEVELS)  ! Tiled soil temperature
      INTEGER :: JSMCL_TILE(SM_LEVELS)   ! Tiled soil moisture content in layers
! Not needed MRD
!!!      INTEGER :: JSTHU_TILE(SM_LEVELS)   ! Tiled unfrozen soil moisture fraction
      INTEGER :: JSTHF_TILE(SM_LEVELS)   ! Tiled frozen soil moisture fraction
      INTEGER :: JSNOW_DEPTH3L(3)        ! Tiled snow depth
      INTEGER :: JSNOW_MASS3L(3)         ! Tiled snow mass
      INTEGER :: JSNOW_TMP3L(3)          ! Tiled snow temperature
      INTEGER :: JSNOW_RHO3L(3)          ! Tiled snow density
      INTEGER :: JSNOW_RHO1L             ! Tiled mean snow density
      INTEGER :: JSNOW_AGE               ! Tiled snow age
      INTEGER :: JSNOW_FLG3L             ! Flag for use of 3 level snow scheme
! Lestevens Sept 2012 - adding new progs for CASA-CNP
      INTEGER :: JCPOOL_TILE(10)
      INTEGER :: JNPOOL_TILE(10)
      INTEGER :: JPPOOL_TILE(12)
      INTEGER :: JSOIL_ORDER
      INTEGER :: JNIDEP
      INTEGER :: JNIFIX
      INTEGER :: JPWEA
      INTEGER :: JPDUST
      INTEGER :: JGLAI
      INTEGER :: JPHENPHASE
      INTEGER :: JCABLE_LAI(1)   

! Addresses in D1 array of primary variables: scalars
      COMMON/CARGPT_ATMOS/                                              &
! Data variables
     &  JTSTAR, JLAND, JPSTAR, JTSTAR_ANOM,                             &
     &  JFRAC_LAND, JTSTAR_LAND, JTSTAR_SEA, JTSTAR_SICE,               &
     &  JSICE_ALB, JLAND_ALB,                                           &
      ! Cloud fields
     &  JCCB, JCCT, JCCLWP,                                             &
      ! Boundary layer fields
     &  JZH, JT1_SD, JQ1_SD, JNTML, JNTDSC, JNBDSC, JCUMULUS,           &
      ! Soil Ancillary fields
     &  JSAT_SOILW_SUCTION, JTHERM_CAP, JTHERM_COND,                    &
     &  JVOL_SMC_CRIT, JVOL_SMC_WILT, JVOL_SMC_SAT,                     &
     &  JSAT_SOIL_COND, JCLAPP_HORN,                                    &
      ! Vegetation Ancillary fields
     &  JCANOPY_WATER, JSURF_CAP, JSURF_RESIST, JROOT_DEPTH,            &
     &  JINFILT, JVEG_FRAC, JLAI, JCANHT,                               &
     &  JZ0, JSFA, JMDSA, JGS,                                          &
      ! Orographic Ancillary fields
     &  JOROG, JOROG_SD, JOROG_SIL, JOROG_HO2,                          &
     &  JOROG_GRAD_X, JOROG_GRAD_Y,                                     &
     &  JOROG_GRAD_XX, JOROG_GRAD_XY, JOROG_GRAD_YY,                    &
      ! Sea/Sea Ice
     &  JU_SEA, JV_SEA, JU_0_P, JV_0_P,                                 &
     &  JICE_FRACTION, JICE_THICKNESS, JTI,                             &
     &  JICE_FRACT_CAT, JICE_THICK_CAT, JTI_CAT,                        &
      ! Snow
     &  JSNODEP, JSNODEP_SEA, JSNSOOT, JCATCH_SNOW, JSNOW_GRND,         &
     &  JSNODEP_SEA_CAT,                                                &

      ! Aerosol emission fields,including mineral dust parent soil props
     &  JSOIL_CLAY, JSOIL_SILT, JSOIL_SAND,JDUST_MREL1,JDUST_MREL2,     &
     &  JDUST_MREL3,JDUST_MREL4,JDUST_MREL5,JDUST_MREL6,                &
     &  JSO2_EM, JDMS_EM, JSO2_HILEM, JNH3_EM, JSOOT_EM, JSOOT_HILEM,   &
     &  JBMASS_EM, JBMASS_HILEM, JOCFF_EM, JOCFF_HILEM, JDMS_CONC,      &
     &  JDMS_OFLUX,                                                     &
      ! User ancillary fields
     &  JUSER_ANC1,  JUSER_ANC2, JUSER_ANC3, JUSER_ANC4,                &
     &  JUSER_ANC5,  JUSER_ANC6, JUSER_ANC7, JUSER_ANC8,                &
     &  JUSER_ANC9,  JUSER_ANC10, JUSER_ANC11, JUSER_ANC12,             &
     &  JUSER_ANC13,  JUSER_ANC14, JUSER_ANC15, JUSER_ANC16,            &
     &  JUSER_ANC17,  JUSER_ANC18, JUSER_ANC19, JUSER_ANC20,            &
      ! Tracer fluxes - kdcorbin, 04/10
     &  JTRACER_FLUX1, JTRACER_FLUX2, JTRACER_FLUX3, JTRACER_FLUX4,     &
     &  JTRACER_FLUX5, JTRACER_FLUX6, JTRACER_FLUX7, JTRACER_FLUX8,     &
     &  JTRACER_FLUX9, JTRACER_FLUX10,JTRACER_FLUX11,JTRACER_FLUX12,    &
     &  JTRACER_FLUX13,JTRACER_FLUX14,JTRACER_FLUX15,JTRACER_FLUX16,    &
     &  JTRACER_FLUX17,JTRACER_FLUX18,JTRACER_FLUX19,JTRACER_FLUX20,    &
     &  JTSLAB,JTCLIM,JHCLIM,JCHEAT,JOIFLX,JUICE,JVICE,                 &
     &  JSIG11NE,JSIG11SE,JSIG11SW,JSIG11NW,                            &
     &  JSIG12NE,JSIG12SE,JSIG12SW,JSIG12NW,                            &
     &  JSIG22NE,JSIG22SE,JSIG22SW,JSIG22NW,                            &
      ! Pointers for ATMOSPHERE model constants. Scalars only.
     &  JETATHETA, JETARHO, JRHCRIT, JSOIL_THICKNESS,                   &
     &  Jzseak_theta,Jck_theta,Jzseak_rho,Jck_rho,                      &
      ! pointers for input variable grid info.
     &  JLAMBDA_INPUT_P,  JLAMBDA_INPUT_U,                              &
     &  JPHI_INPUT_P, JPHI_INPUT_V,                                     &
      ! pointers for tiled vegetation and triffid
     &  JFRAC_TYP, JLAI_PFT, JCANHT_PFT, JDISTURB, JSOIL_ALB,           &
     &  JFRAC_CON1, JFRAC_CON2, JFRAC_CON3, JFRAC_CON4, JFRAC_CON5,     &
     &  JFRAC_CON6, JFRAC_CON7, JFRAC_CON8, JFRAC_CON9,                 &
     &  JSOIL_CARB,                                                     &
     &  JSOIL_CARB1, JSOIL_CARB2, JSOIL_CARB3, JSOIL_CARB4,             &
     &  JNPP_PFT_ACC, JG_LF_PFT_ACC, JG_PHLF_PFT_ACC,                   &
     &  JRSP_W_PFT_ACC, JRSP_S_ACC,                                     &
     &  JRSP_S_ACC1, JRSP_S_ACC2, JRSP_S_ACC3,                          &
     &  JRSP_S_ACC4, JCAN_WATER_TILE, JCATCH_TILE,                      &
     &  JINFIL_TILE, JRGRAIN_TILE, JSNODEP_TILE, JTSTAR_TILE, JZ0_TILE, &
     &  JDOLR, JLW_DOWN, JSW_TILE,JNET_FLUX,JNET_MFLUX,                 &
      ! pointers for carbon cycle
     &  J_CO2FLUX,J_CO2_EMITS,                                          &
      ! pointers for large-scale hydrology
     &  JFEXP, JTI_MEAN, JTI_SIG, JGAMMA_INT,                           &
     &  JWATER_TABLE, JFSFC_SAT, JF_WETLAND, JSTHZW,                    &
     &  JA_FSAT,      JC_FSAT,   JA_FWET,    JC_FWET,                   &
      ! pointers for river routing
     &  JRIV_SEQUENCE, JRIV_DIRECTION, JRIV_STORAGE,                    &
     &  JTOT_SURFROFF, JTOT_SUBROFF, JRIV_INLANDATM                     &
      ! pointers for coupling fields
     & , JC_SOLAR, JC_BLUE, JC_DOWN, JC_LONGWAVE, JC_TAUX               &
     & , JC_TAUY, JC_WINDMIX, JC_SENSIBLE, JC_SUBLIM                    &
     & , JC_EVAP, JC_BOTMELTN, JC_TOPMELTN, JC_LSRAIN                   & 
     & , JC_LSSNOW, JC_CVRAIN, JC_CVSNOW, JC_RIVEROUT,JC_PRESS

! TYPPTRA end
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

! CAOPTR Holds address pointers for atmosphere-to-ocean coupling fields
! required by SWAP_A2O and SWAP_O2A.
! 5.5 28/02/03 Add pointers for river outflow on atmos grid
! later remove runoff and ocenpts pointers  C.Bunton
! 6.2 24/02/06 Add pointers for 10m windspeed            J.Gunson
! 6.2 29/11/05 Add pointers for ice velocities on OCEAN grid  A.B.Keen
!
! 6.2 24/02/06 Add pointers for ocean DMS flux to atmosphere    J.Gunson
! Pointers needed by SWAP_A2O and SWAP_O2A

      INTEGER :: JO_SST    ! Sea-surface temperature on ocean grid
      INTEGER :: JO_UCURR  ! Surface zonal current on ocean grid
      INTEGER :: JA_AICE   ! Fractional ice conc. on atmos grid
      INTEGER :: JA_AICEN  ! Category frac ice conc. on atmos grid
      INTEGER :: JO_AICEN  ! Category frac ice conc. on ocean grid
      COMMON /AOPTR/JO_SST,JO_UCURR,JA_AICE,JA_AICEN,JO_AICEN

      ! Pointers needed by SWAP_A2O

      INTEGER:: JA_TAUX       ! Surface x-windstress on atmos grid
      INTEGER:: JO_TAUX       ! Surface x-windstress on ocean grid
      INTEGER:: JA_TAUY       ! Surface y-windstress on atmos grid
      INTEGER:: JO_TAUY       ! Surface y-windstress on ocean grid
      INTEGER:: JA_WINDMIX    ! Windmixing power on atmos grid
      INTEGER:: JO_WINDMIX    ! Windmixing power on ocean grid
      INTEGER:: JA_W10        ! 10m windspeed
      INTEGER:: JO_W10
      INTEGER:: JA_DDDL1_D1   ! Dust dry dep lyr 1
      INTEGER:: JA_DDDL1_D2
      INTEGER:: JA_DDDL1_D3
      INTEGER:: JA_DDDL1_D4
      INTEGER:: JA_DDDL1_D5
      INTEGER:: JA_DDDL1_D6
      INTEGER:: JA_DDDL2_D1   ! Dust dry dep lyr 2
      INTEGER:: JA_DDDL2_D2
      INTEGER:: JA_DDDL2_D3
      INTEGER:: JA_DDDL2_D4
      INTEGER:: JA_DDDL2_D5
      INTEGER:: JA_DDDL2_D6
      INTEGER:: JA_DWDLS_D1   ! Dust wet dep lrg scl rain
      INTEGER:: JA_DWDLS_D2
      INTEGER:: JA_DWDLS_D3
      INTEGER:: JA_DWDLS_D4
      INTEGER:: JA_DWDLS_D5
      INTEGER:: JA_DWDLS_D6
      INTEGER:: JA_DWDCN_D1   ! Dust wet dep convctv rain
      INTEGER:: JA_DWDCN_D2
      INTEGER:: JA_DWDCN_D3
      INTEGER:: JA_DWDCN_D4
      INTEGER:: JA_DWDCN_D5
      INTEGER:: JA_DWDCN_D6
      INTEGER:: JO_DDEPD1     ! Dust deposition : ocean
      INTEGER:: JO_DDEPD2
      INTEGER:: JO_DDEPD3
      INTEGER:: JO_DDEPD4
      INTEGER:: JO_DDEPD5
      INTEGER:: JO_DDEPD6
      INTEGER:: JA_SOLAR      ! Net downward SW at surf on atmos grid
      INTEGER:: JA_BLUE       ! Net blueband SW at surf on atmos grid
      INTEGER:: JO_BLUE       ! Net blueband SW at surf on ocean grid
      INTEGER:: JA_EVAP       ! Net evaporation over sea on atmos grid
      INTEGER:: JA_LONGWAVE   ! Net downward LW at surf on atmos grid
      INTEGER:: JA_SENSIBLE  ! Sensible heat flux over sea on atmos grid
      INTEGER:: JO_HEATFLUX   ! Non penetrative heatflux on ocean grid
      INTEGER:: JA_LSSNOW     ! Large-scale snowfall rate on atmos grid
      INTEGER:: JA_CVSNOW     ! Convective snowfall rate on atmos grid
      INTEGER:: JA_LSRAIN     ! Large-scale rainfall rate on atmos grid
      INTEGER:: JA_CVRAIN     ! Convective rainfall rate on atmos grid
      INTEGER:: JO_PMINUSE    ! Precipitation-evaporation on ocean grid
      INTEGER:: JA_SLOWRUNOFF ! Slow (sub-surface) runoff on atmos grid
      INTEGER:: JA_FASTRUNOFF ! Fast (surface) runoff on atmos grid
      INTEGER:: JA_RIVEROUT   ! Total river outflow on atmos grid
      INTEGER:: JA_OCENTPTS   ! Ocean entry point index to atmos landpts
      INTEGER:: JO_RIVEROUT   ! Total river outflow on ocean grid
      INTEGER:: JA_co2        ! atmos level 1 co2 conc.
      INTEGER:: JO_co2
      INTEGER:: JA_co2flux    ! ocean co2 flux.
      INTEGER:: JO_co2flux
      INTEGER:: JA_dmsflux    ! ocean DMS flux.
      INTEGER:: JO_dmsflux
      INTEGER:: JO_SNOWFALL   ! Snowfall rate on ocean grid
      INTEGER:: JA_SUBLIM     ! Sublimation on atmos grid
      INTEGER:: JO_SUBLIM     ! Sublimation on ocean grid
      INTEGER:: JA_BOTMELT    ! Diffusive heat thro ice on atmos grid
      INTEGER:: JA_BOTMELTN   ! Diff heat thro ice (ncat)
      INTEGER:: JO_BOTMELT    ! Diffusive heat thro ice on ocean grid
      INTEGER:: JA_TOPMELT    ! Seaice top melting flux on atmos grid
      INTEGER:: JA_TOPMELTN   ! Seaice top melting flux (ncat)
      INTEGER:: JO_TOPMELT    ! Seaice top melting flux on ocean grid
      INTEGER:: JO_AICE       ! Sea ice concentration on ocean grid
      INTEGER:: JA_PRESS      ! Surface pressure on atmos grid
      INTEGER :: JC_U10
      INTEGER :: JC_V10
      INTEGER :: JC_CO2
 
      COMMON /A2OPTR/                                                   &
     &  JA_TAUX,JO_TAUX,JA_TAUY,JO_TAUY,JA_WINDMIX,JO_WINDMIX,          &
     &  JA_W10, JO_W10,                                                 &
     &  JA_DDDL1_D1,JA_DDDL1_D2,JA_DDDL1_D3,JA_DDDL1_D4,JA_DDDL1_D5,    &
     &  JA_DDDL1_D6,JA_DDDL2_D1,JA_DDDL2_D2,JA_DDDL2_D3,JA_DDDL2_D4,    &
     &  JA_DDDL2_D5,JA_DDDL2_D6,JA_DWDLS_D1,JA_DWDLS_D2,JA_DWDLS_D3,    &
     &  JA_DWDLS_D4,JA_DWDLS_D5,JA_DWDLS_D6,JA_DWDCN_D1,JA_DWDCN_D2,    &
     &  JA_DWDCN_D3,JA_DWDCN_D4,JA_DWDCN_D5,JA_DWDCN_D6,                &
     &  JO_DDEPD1,JO_DDEPD2,JO_DDEPD3,JO_DDEPD4,JO_DDEPD5,JO_DDEPD6,    &
     &  JA_SOLAR,JA_BLUE,JO_BLUE,JA_EVAP,JA_LONGWAVE,JA_SENSIBLE,       &
     &  JO_HEATFLUX,JA_LSSNOW,JA_CVSNOW,JA_LSRAIN,JA_CVRAIN,JO_PMINUSE, &
     &  JA_SLOWRUNOFF,JA_FASTRUNOFF,JA_OCENTPTS,JA_RIVEROUT,            &
     &  JO_RIVEROUT,JA_co2, JO_co2, JA_co2flux, JO_co2flux,             &
     &  JA_dmsflux, JO_dmsflux,                                         &
     &  JO_SNOWFALL,JA_SUBLIM,JO_SUBLIM,JA_BOTMELT,JO_BOTMELT,          &
     &  JA_TOPMELT,JO_TOPMELT,JA_BOTMELTN,JA_TOPMELTN,JO_AICE,JA_PRESS, &
     &  JC_U10, JC_V10, JC_CO2 

      INTEGER:: JO_TSTAR           ! Surface temperature on ocean grid
      INTEGER:: JA_TSTAR           ! Surface temperature on atmos grid
      INTEGER:: JA_UCURR           ! Surface zonal current on atmos grid
      INTEGER:: JO_VCURR           ! Surface merid current on ocean grid
      INTEGER:: JA_VCURR           ! Surface merid current on atmos grid
      INTEGER:: JO_UICE            ! Zonal seaice vel on ocean grid
      INTEGER:: JO_VICE            ! Merid seaice vel on ocean grid
      INTEGER:: JO_ICEDEPTH        ! Ice depth on ocean grid
      INTEGER:: JA_ICEDEPTH        ! Ice depth on atmos grid
      INTEGER:: JA_ICEDEPTHN       ! Ice depth on atmos grid (ncat)
      INTEGER:: JO_SNOWDEPTH       ! Snow depth on ocean grid
      INTEGER:: JA_SNOWDEPTH       ! Snow depth on atmos grid
      INTEGER:: JA_SNOWDEPTHN      ! Snow depth on atmos grid (ncat)

      COMMON /O2APTR/                                                   &
     &  JO_TSTAR,JA_TSTAR,JA_UCURR,JO_VCURR,JA_VCURR,                   &
     &  JO_UICE,JO_VICE,                                                &
     &  JO_ICEDEPTH,JA_ICEDEPTH,JO_SNOWDEPTH,JA_SNOWDEPTH,              &
     &  JA_ICEDEPTHN,JA_SNOWDEPTHN
! CAOPTR end
! C_KAPPAI start

! Thermal conductivity of sea-ice (W per m per K).
        REAL,PARAMETER:: KAPPAI=2.09

! Thermal conductivity of sea water (W per m per K).
        REAL,PARAMETER:: kappas=0.31

! Snow density (Kg per m**3)
        REAL,PARAMETER:: rhosnow=330.0

! Effective thickness of sea-ice surface layer (m).
        REAL,PARAMETER:: DE = 0.1

! C_KAPPAI end
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
! cntlall.h required for ltimer
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
!*L --------------------- Comdeck: CHISTORY ----------------------------
!LL
!LL  Purpose: COMMON block for history data needed by top level (C0)
!LL           routines, and passed from run to run.  Mostly set by
!LL           the User Interface.
!LL
!LL           Note that CHISTORY *CALLs ALL individual history comdecks
!LL
!LL  Author : A. Sangster
!LL
!LL  Model            Modification history
!LL version  Date
!LL  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!LL                 contents.  RTHBarnes.
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LLEND----------------------------------------------------------------
!*
!CC   *CALL CHSUNITS
! ----------------------- Comdeck: IHISTO   ----------------------------
! Description: COMDECK defining Integer History variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!

      ! Array containing model data time (Same as MODEL_BASIS_TIME/MODEL
      ! ANALYSIS_HRS depending whether before/after assimilation)
      INTEGER :: MODEL_DATA_TIME(6)

      ! Indicator for next mean period to be processed
      INTEGER :: RUN_MEANCTL_RESTART

      ! Indicator of operational run type
      INTEGER :: RUN_INDIC_OP

      ! Final target date for the run
      INTEGER :: RUN_RESUBMIT_TARGET(6)

      ! Last field written/read per FT unit
      INTEGER ::FT_LASTFIELD(20:NUNITS)

! History Common Block for overall model integers variables.

      COMMON /IHISTO/                                                   &
     &  MODEL_DATA_TIME,                                                &
     &  RUN_MEANCTL_RESTART, RUN_INDIC_OP,                              &
     &  RUN_RESUBMIT_TARGET, FT_LASTFIELD

      NAMELIST /NLIHISTO/                                               &
     &  MODEL_DATA_TIME,                                                &
     &  RUN_MEANCTL_RESTART, RUN_INDIC_OP,                              &
     &  RUN_RESUBMIT_TARGET, FT_LASTFIELD

! IHISTO end
! ----------------------- Comdeck: CHISTO   ----------------------------
! Description: COMDECK defining Character History variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.1  18/04/96  Add RUN_IN for qxhistreport.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
      CHARACTER(LEN=10) :: RUN_HIST_TYPE       ! Type of history file
      CHARACTER*8  RUN_TYPE            ! Type of run
      CHARACTER*14 RUN_COMPCODE        ! Run completion code
      CHARACTER*14 RUN_LAST_MEAN       ! Last mean dump created by run
      ! Appears Unused                 ! for pp fields
      CHARACTER*1  RUN_MEANS_TO_DO     ! Flag indicating the run stopped
                                       ! before creating next mean dump
      CHARACTER*1  RUN_OCEAN_FIRST     ! Flag set to true if ocean to be
                                       ! run first
      CHARACTER*8  RUN_JOB_NAME        ! Jobname this run
      CHARACTER*5  RUN_ID              ! Expt./Job id for this run
      CHARACTER*1  RUN_RESUBMIT        ! Flag controlling auto resubmit
      CHARACTER*12 RUN_RESUBMIT_Q      ! Job queue to which resubmit run
      CHARACTER*20 RUN_RESUBMIT_TIME   ! Time at which run resubmits
      CHARACTER*6  RUN_RESUBMIT_CPU    ! Time limit for resubmitted job
      CHARACTER*6  RUN_RESUBMIT_MEMORY ! Resubmitted job's memory limit
      CHARACTER*2  RUN_RESUBMIT_PRTY   ! Resubmitted job intra q prty
      CHARACTER*8  RUN_RESUBMIT_JOBNAME! Resubmitted jobname
      CHARACTER*1  FT_ACTIVE(20:NUNITS) ! "Y" if file partly written

      ! History Common Block for overall model character variables.

      COMMON /CHISTO/                                                   &
     &  RUN_HIST_TYPE, RUN_TYPE, RUN_COMPCODE, RUN_LAST_MEAN,           &
     &  RUN_MEANS_TO_DO, RUN_OCEAN_FIRST, RUN_JOB_NAME, RUN_ID,         &
     &  RUN_RESUBMIT, RUN_RESUBMIT_Q, RUN_RESUBMIT_TIME,                &
     &  RUN_RESUBMIT_CPU, RUN_RESUBMIT_MEMORY, RUN_RESUBMIT_PRTY,       &
     & RUN_RESUBMIT_JOBNAME, FT_ACTIVE

      NAMELIST /NLCHISTO/                                               &
     &  RUN_HIST_TYPE, RUN_TYPE, RUN_COMPCODE, RUN_LAST_MEAN,           &
     &  RUN_MEANS_TO_DO, RUN_OCEAN_FIRST, RUN_JOB_NAME, RUN_ID,         &
     &  RUN_RESUBMIT, RUN_RESUBMIT_Q, RUN_RESUBMIT_TIME,                &
     &  RUN_RESUBMIT_CPU, RUN_RESUBMIT_MEMORY, RUN_RESUBMIT_PRTY,       &
     & RUN_RESUBMIT_JOBNAME, FT_ACTIVE

! CHISTO end
! ----------------------- Comdeck: IHISTG   ----------------------------
! Description: COMDECK defining Integer History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      ! No. of tsteps completed this run
      INTEGER :: LENGTH(N_INTERNAL_MODEL_MAX)

      ! Model end time this run
      INTEGER :: ACTUAL_ENDT(6,N_INTERNAL_MODEL_MAX)

      ! These 2 appears to be purely diagnostic, and not really used.

      ! History block copy of A/O_STEP held in file CTIME
      INTEGER :: H_STEPim(N_INTERNAL_MODEL_MAX)

      ! No of steps in coupling period
      INTEGER :: H_GROUPim(N_INTERNAL_MODEL_MAX)

      ! No of means activated
      INTEGER :: MEAN_OFFSETim(N_INTERNAL_MODEL_MAX)

      ! Offset between MEAN_REFTIME and model basis time(in model dumps)
      INTEGER :: OFFSET_DUMPSim(N_INTERNAL_MODEL_MAX)

      ! No of mean periods chosen
      INTEGER :: MEAN_NUMBERim(N_INTERNAL_MODEL_MAX)

      ! Indicators used to correct logical units are used for
      ! atmos/ocean partial sum dump I/O
      INTEGER :: RUN_MEANCTL_INDICim(4,N_INTERNAL_MODEL_MAX)

      ! History Common Block for generic model integer variables.

      COMMON /IHISTG/                                                   &
     &  H_STEPim, H_GROUPim, MEAN_OFFSETim, OFFSET_DUMPSim,             &
     & MEAN_NUMBERim, RUN_MEANCTL_INDICim

      NAMELIST /NLIHISTG/                                               &
     &  H_STEPim, H_GROUPim, MEAN_OFFSETim, OFFSET_DUMPSim,             &
     & MEAN_NUMBERim, RUN_MEANCTL_INDICim

! IHISTG end
! ----------------------- Comdeck: CHISTG   ----------------------------
! Description: COMDECK defining Character History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.4  30/05/97  Added vars LASTATMim, CURRATMim, LASTDMPim.  K Rogers
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      CHARACTER*14 END_DUMPim(N_INTERNAL_MODEL_MAX)!most recent dumpname
      CHARACTER*80 RESTARTim(N_INTERNAL_MODEL_MAX) !current restart dump
      CHARACTER*14 SAFEDMPim(N_INTERNAL_MODEL_MAX)
! Name of old safe restart dump
      CHARACTER*14 NEWSAFEim(N_INTERNAL_MODEL_MAX)
! Name of new safe restart dump
      CHARACTER*14 LASTATMim(N_INTERNAL_MODEL_MAX) ! Keep name of last
!                                                  ! atmos restart dump
!                                                  ! until ocean dump
      CHARACTER*14 CURRATMim(N_INTERNAL_MODEL_MAX) ! Keep name of
!                                                  ! current atmos
!                                                  ! restart dump
      CHARACTER*14 LASTDMPim(N_INTERNAL_MODEL_MAX) ! Keep name of last
!                                                  ! atmos/ocean dumps
!                                                  ! until meaning done

!
!
! History Common Block for generic model characters variables.
!
      COMMON /CHISTG/                                                   &
     &  END_DUMPim, RESTARTim,                                          &
     &  SAFEDMPim, NEWSAFEim, LASTATMim, CURRATMim, LASTDMPim

      NAMELIST /NLCHISTG/                                               &
     &  END_DUMPim, RESTARTim,                                          &
     &  SAFEDMPim, NEWSAFEim, LASTATMim, CURRATMim, LASTDMPim

! CHISTG end
!*L --------------------- Comdeck: CLFHIST  ----------------------------
!LL
!LL  Purpose: COMDECK defining unit numbers relevant to history file
!LL           and variables used to hold the logical to physical
!LL           file associations made within the model
!LL
!LL  Author : A. Sangster
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LL                  Version 5  18/6/90
!LL
!LL  Model             Modification history from model version 3.0
!LL version  Date
!LL
!LL  3.4  30/09/94  Add files MURKFILE,OUSRANCL,OUSRMULT at 109,113,114
!LL  3.4  05/09/94  Add files USRANCIL,USRMULTI at unit nos. 111,112.
!LL
!LL  3.3  22/11/93  Add file SOURCES at unit number 110. R.T.H.Barnes.
!LL 3.2     28/05/93  Add file BAS_IND at unit number 58. M.Carter.
!LL  Vn3.0  12/02/93 - Variables PERTURB and TRANSP equivalenced to unit
!LL                    numbers 37, and 97 respectively. C.S. Douglas
!LL  3.4  1/8/94     Revised Obs file specification: Stuart Bell
!LL  3.5  01/05/95  Sub-models stage 1: History/control files. RTHBarnes
!    4.0  22/09/95  Added units for Spectral data for Radiation scheme.
!                                        (J. M. Edwards)
!LL  4.1  11/03/96  Introduce Wave sub-model.  RTHBarnes.
!    4.1  26/02/96  Associate new env. variables SO2NATEM and CHEMOXID
!                   with unit nos. 115 & 116. Rename SOURCES to
!                   SULPEMIS. D. Robinson.
!  4.3   18/3/97  Add aerosol forcings of climate change.  Will Ingram
!  4.4   4/7/97   Add ANLINCR  Chris Jones/Stuart Bell
!LL  4.4   12/9/97  Associate ancillary file EVs for initial surface
!LL                 type fracs, initial vegetation state and vegetation
!LL                 disturbance with unit no.s 135-137 R. Betts
!LL  4.4  17/10/97  Associate env var. CACHED with Unit 138. D Robinson
!LL  4.5  22/04/98  Add new ancillary file for soot emissions:
!LL                 SOOTEMIS - in I/O unit 139. R.Rawlins
!LL  4.5  29/07/98  Add new variables ALABCOU5/6/7/8. D. Robinson.
!LL  4.5  17/08/98  Add new variables OLABCOU1/2/3/4. Remove
!LL                 OLABCOUT. D. Robinson.
!LL  5.1  13/04/00  TDF_dump added. ANLINCR changed to IAU_inc.
!LL                 Adam Clayton
!LL  5.2  21/08/00  Add an extra op macro for VAR plus 1 user pp
!LL                 output stream. R Rawlins
!LL  5.2  18/01/01  Add VERT_LEV and attach to unit 90. D. Robinson
!LL  5.3  26/10/01  Add LANDFRAC and attach to unit 120.
!LL                 Free units 75-79 (OBS06-OBS10). D. Robinson
!    5.3  24/10/01  Add IDEALISE and attach to unit 106. A. Malcolm
!    5.3  14/11/00  Added TPPSOZON (tropopause-based ozone). Dave Tan
!    5.4  29/08/02  Add ALABCIN1/2 & attach to unit 125/126. D Robinson
!    5.5  17/02/03  Add Wave model boundary & interface files
!                                                 D.Holmes-Bell
!    5.5  30/12/02  Add DUSTSOIL/BIOMASS and attach to unit 75/76.
!                   RIVSTOR/RIVSEQ/RIVDIR to units 77-79.
!                   D Robinson
!    6.1  07/04/04  Add DMSCONC to unit 95.        A. Jones
!    6.1  08/11/04  Alter names of River Routing files. R.Sharp
!    6.2  26/01/06  Include iteration count file name string
!                                                    T. Edwards
!    6.2 24/00/05   Change STRATOUT to PPSMC and MESOUT to PPSCREEN,
!                   files for soil moisture nudging scheme. Clive Jones
!    6.2  17/03/06  Add SURFEMIS/AIRCREMS/STRATEMS/EXTRAEMS/RADONEMS
!                   filenames to unit nos 130-134 for UKCA emissions
!                   files. WINITIAL/WSTART/WRESTART/WAVANL/WAVANCIN
!                   removed. F. O'Connor
!LL
!LL  Type declarations
!LL
!LL
!LL  Logical Filenames used in the model
!LL
      CHARACTER*80 HKFILE,PPXREF,CONFIG,STASHCTL,NAMELIST,OUTPUT,       &
     &             OUTPUT2,MCTL,ICTL,PHIST,IHIST,THIST,FTXX,            &
     &             CACHE1,CACHE2,ASWAP,OSWAP,AOTRANS,                   &
     &             AINITIAL,ASTART,ARESTART,AOPSUM1,AOPSUM2,AOPSUM3,    &
     &             AOPSUM4,AOMEAN,SSU,                                  &
     &             OZONE,SMCSNOWD,DSOILTMP,SOILTYPE,VEGTYPE,SSTIN,      &
     &             SICEIN,PERTURB,MASK,                                 &
     &             OINITIAL,OSTART,ORESTART,AOPSTMP1,AOPSTMP2,AOPSTMP3, &
     &             AOPSTMP4,                                            &
     &             WFIN,HFLUXIN,PMEIN,ICEFIN,AIRTMP,                    &
     &             SWSPECTD,                                            &
     &             PP0,PP1,PP2,PP3,PP4,PP5,PP6,PP7,PP8,PP9,             &
     &             PPVAR,PP10,                                          &
     &             OBS01,OBS02,OBS03,OBS04,OBS05,                       &
     &             DUSTSOIL,BIOMASS,RIVSTOR,RIVCHAN,RIVER2A,            &
     &             SURFEMIS, AIRCREMS, STRATEMS, EXTRAEMS, RADONEMS,    &
     &             LWSPECTD,WAVEOUT,SURGEOUT,PPSCREEN,PPSMC,WFOUT,      &
     &          HFLUXOUT,FLXCROUT,PMEOUT,ICEFOUT,MOSOUT,SSTOUT,SICEOUT, &
     &             CURNTOUT,DMSCONC,OROG,OLABCIN,OCNDEPTH,CURNTIN,      &
     &             FLUXCORR,SLABHCON,ATMANL,OCNANL,BAS_IND,             &
     &             TRANSP,ATRACER,OTRACER,SULPEMIS,USRANCIL,USRMULTI,   &
     &             OUSRANCL,OUSRMULT,MURKFILE,                          &
     &             ALABCIN1,ALABCIN2,                                   &
     &             ALABCOU1,ALABCOU2,ALABCOU3,ALABCOU4,                 &
     &             ALABCOU5,ALABCOU6,ALABCOU7,ALABCOU8,CARIOLO3,        &
     &             OLABCOU1,OLABCOU2,OLABCOU3,OLABCOU4,                 &
     &             WLABCOU1,WLABCOU2,WLABCOU3,WLABCOU4,HORZGRID,        &
     &             TDF_dump,IAU_inc,                                    &
     &             LANDFRAC,                                            &
     &             SO2NATEM,CHEMOXID,AEROFCG,FRACINIT,VEGINIT,DISTURB,  &
     &             CACHED,SOOTEMIS,                                     &
     &             CO2EMITS,TPPSOZON,                                   &
     &             VERT_LEV,VAR_GRID,                                   &
     &             IDEALISE,ICFILE,                                     &
     &             ARCLBIOG,ARCLBIOM,ARCLBLCK,ARCLSSLT,ARCLSULP,        &
     &             ARCLDUST,ARCLOCFF,ARCLDLTA,RPSEED,OCFFEMIS

!
      CHARACTER*80 MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                               ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        PHIST_UNIT,                                               &
                                 ! Permanent history file unit
     &        IHIST_UNIT,                                               &
                                 ! Interim history file unit
     &        THIST_UNIT,                                               &
                                 ! Temporary history file unit
     &        FTXX_UNIT,                                                &
                                 ! Logical/physical file associations
     &        HKFILE_UNIT        ! Operational houskeeping file unit
!*
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(PHIST_UNIT =10)
      PARAMETER(IHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)
      PARAMETER(FTXX_UNIT  =13)
!
! Namelist of all permissible logical files.
!
      NAMELIST / NLCFILES /                                             &
     &             HKFILE,PPXREF,CONFIG,STASHCTL,NAMELIST,OUTPUT,       &
     &             OUTPUT2,MCTL,ICTL,PHIST,IHIST,THIST,FTXX,            &
     &             CACHE1,CACHE2,ASWAP,OSWAP,AOTRANS,                   &
     &             AINITIAL,ASTART,ARESTART,AOPSUM1,AOPSUM2,AOPSUM3,    &
     &             AOPSUM4,AOMEAN,SSU,                                  &
     &             OZONE,SMCSNOWD,DSOILTMP,SOILTYPE,VEGTYPE,SSTIN,      &
     &             SICEIN,PERTURB,MASK,                                 &
     &             OINITIAL,OSTART,ORESTART,AOPSTMP1,AOPSTMP2,AOPSTMP3, &
     &             AOPSTMP4,                                            &
     &             WFIN,HFLUXIN,PMEIN,ICEFIN,AIRTMP,                    &
     &             SWSPECTD,                                            &
     &             PP0,PP1,PP2,PP3,PP4,PP5,PP6,PP7,PP8,PP9,             &
     &             PPVAR,PP10,                                          &
     &             OBS01,OBS02,OBS03,OBS04,OBS05,                       &
     &             DUSTSOIL,BIOMASS,RIVSTOR,RIVCHAN,RIVER2A,            &
     &             SURFEMIS, AIRCREMS, STRATEMS, EXTRAEMS, RADONEMS,    &
     &             LWSPECTD,WAVEOUT,SURGEOUT,PPSCREEN,PPSMC,WFOUT,      &
     &          HFLUXOUT,FLXCROUT,PMEOUT,ICEFOUT,MOSOUT,SSTOUT,SICEOUT, &
     &             CURNTOUT,DMSCONC,OROG,OLABCIN,OCNDEPTH,CURNTIN,      &
     &             FLUXCORR,SLABHCON,ATMANL,OCNANL,BAS_IND,             &
     &             TRANSP,ATRACER,OTRACER,SULPEMIS,USRANCIL,USRMULTI,   &
     &             OUSRANCL,OUSRMULT,MURKFILE,                          &
     &             ALABCIN1,ALABCIN2,                                   &
     &             ALABCOU1,ALABCOU2,ALABCOU3,ALABCOU4,                 &
     &             ALABCOU5,ALABCOU6,ALABCOU7,ALABCOU8,CARIOLO3,        &
     &             OLABCOU1,OLABCOU2,OLABCOU3,OLABCOU4,                 &
     &             WLABCOU1,WLABCOU2,WLABCOU3,WLABCOU4,HORZGRID,        &
     &             TDF_dump,IAU_inc,                                    &
     &             LANDFRAC,                                            &
     &             SO2NATEM,CHEMOXID,AEROFCG,FRACINIT,VEGINIT,DISTURB,  &
     &             CACHED,SOOTEMIS,                                     &
     &             CO2EMITS,TPPSOZON,                                   &
     &             VERT_LEV,VAR_GRID,                                   &
     &             IDEALISE,ICFILE,                                     &
     &             ARCLBIOG,ARCLBIOM,ARCLBLCK,ARCLSSLT,ARCLSULP,        &
     &             ARCLDUST,ARCLOCFF,ARCLDLTA,RPSEED,OCFFEMIS

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(PHIST      ,MODEL_FT_UNIT(10) ), &
     &(IHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(FTXX      ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &(AOTRANS   ,MODEL_FT_UNIT(17) ),(ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(VEGTYPE   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
     &                                (LWSPECTD   ,MODEL_FT_UNIT(80) ), &
     &(WAVEOUT   ,MODEL_FT_UNIT(81) ),(SURGEOUT   ,MODEL_FT_UNIT(82) ), &
     &(PPSCREEN  ,MODEL_FT_UNIT(83) ),(PPSMC      ,MODEL_FT_UNIT(84) ), &
     &(WFOUT     ,MODEL_FT_UNIT(85) ),(HFLUXOUT   ,MODEL_FT_UNIT(86) ), &
     &(PMEOUT    ,MODEL_FT_UNIT(87) ),(ICEFOUT    ,MODEL_FT_UNIT(88) ), &
     &(MOSOUT    ,MODEL_FT_UNIT(89) ),(VERT_LEV   ,MODEL_FT_UNIT(90) ), &
     &(SSTOUT    ,MODEL_FT_UNIT(91) ),(SICEOUT    ,MODEL_FT_UNIT(92) ), &
     &(CURNTOUT  ,MODEL_FT_UNIT(93) ),(FLXCROUT   ,MODEL_FT_UNIT(94) ), &
     &(DMSCONC   ,MODEL_FT_UNIT(95) ),(OROG       ,MODEL_FT_UNIT(96) ), &
     &(TRANSP    ,MODEL_FT_UNIT(97) ),(OLABCIN    ,MODEL_FT_UNIT(98) ), &
     &(OCNDEPTH  ,MODEL_FT_UNIT(99) ),                                  &
     &(OLABCOU1  ,MODEL_FT_UNIT(100)),(OLABCOU2   ,MODEL_FT_UNIT(101)), &
     &(OLABCOU3  ,MODEL_FT_UNIT(102)),(OLABCOU4   ,MODEL_FT_UNIT(103)), &
     &(IDEALISE  ,MODEL_FT_UNIT(106)),(TDF_dump   ,MODEL_FT_UNIT(107)), &
     &(IAU_inc   ,MODEL_FT_UNIT(108)),(MURKFILE   ,MODEL_FT_UNIT(109)), &
     &(SULPEMIS  ,MODEL_FT_UNIT(110)),(USRANCIL   ,MODEL_FT_UNIT(111)), &
     &(USRMULTI  ,MODEL_FT_UNIT(112)),(OUSRANCL   ,MODEL_FT_UNIT(113)), &
     &(OUSRMULT  ,MODEL_FT_UNIT(114)),(SO2NATEM   ,MODEL_FT_UNIT(115)), &
     &(CHEMOXID  ,MODEL_FT_UNIT(116)),(AEROFCG    ,MODEL_FT_UNIT(117)), &
     &(CO2EMITS  ,MODEL_FT_UNIT(118)),(TPPSOZON   ,MODEL_FT_UNIT(119)), &
     &(LANDFRAC  ,MODEL_FT_UNIT(120)),(WLABCOU1   ,MODEL_FT_UNIT(121)), &
     &(WLABCOU2  ,MODEL_FT_UNIT(122)),(WLABCOU3   ,MODEL_FT_UNIT(123)), &
     &(WLABCOU4  ,MODEL_FT_UNIT(124)),(ALABCIN1   ,MODEL_FT_UNIT(125)), &
     &(ALABCIN2  ,MODEL_FT_UNIT(126)),                                  &
     &(OCFFEMIS  ,MODEL_FT_UNIT(128)),(HORZGRID   ,MODEL_FT_UNIT(129)), &
     &(SURFEMIS  ,MODEL_FT_UNIT(130)),(AIRCREMS   ,MODEL_FT_UNIT(131)), &
     &(STRATEMS  ,MODEL_FT_UNIT(132)),(EXTRAEMS   ,MODEL_FT_UNIT(133)), &
     &(RADONEMS  ,MODEL_FT_UNIT(134)),(FRACINIT   ,MODEL_FT_UNIT(135)), &
     &(VEGINIT   ,MODEL_FT_UNIT(136)),(DISTURB    ,MODEL_FT_UNIT(137)), &
     &(CACHED    ,MODEL_FT_UNIT(138)),(SOOTEMIS   ,MODEL_FT_UNIT(139)), &
     &(ALABCOU1  ,MODEL_FT_UNIT(140)),(ALABCOU2   ,MODEL_FT_UNIT(141)), &
     &(ALABCOU3  ,MODEL_FT_UNIT(142)),(ALABCOU4   ,MODEL_FT_UNIT(143)), &
     &(ALABCOU5  ,MODEL_FT_UNIT(144)),(ALABCOU6   ,MODEL_FT_UNIT(145)), &
     &(ALABCOU7  ,MODEL_FT_UNIT(146)),(ALABCOU8   ,MODEL_FT_UNIT(147)), &
     &(CARIOLO3  ,MODEL_FT_UNIT(148)),(RPSEED     ,MODEL_FT_UNIT(149)), &
     &(PPVAR     ,MODEL_FT_UNIT(150)),(PP10       ,MODEL_FT_UNIT(151)), &
     &(ICFILE    ,MODEL_FT_UNIT(152)),(VAR_GRID   ,MODEL_FT_UNIT(153)), &
     &(ARCLBIOG  ,MODEL_FT_UNIT(154)),(ARCLBIOM   ,MODEL_FT_UNIT(155)), &
     &(ARCLBLCK  ,MODEL_FT_UNIT(156)),(ARCLSSLT   ,MODEL_FT_UNIT(157)), &
     &(ARCLSULP  ,MODEL_FT_UNIT(158)),(ARCLDUST   ,MODEL_FT_UNIT(159)), &
     &(ARCLOCFF  ,MODEL_FT_UNIT(160)),(ARCLDLTA   ,MODEL_FT_UNIT(161))
!
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
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
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

  !     Subroutine arguments
  LOGICAL :: GET_STEP       ! Whether to transfer incoming
  ! data to UM atmos D1 or just ignore it.
  CHARACTER*(*) :: cmessage ! OUT - Error return message

  !     Local variables
  INTEGER :: i, j, k, n           ! Loop counters
  INTEGER :: r_cols, r_rows ! Receive buffer dimensions

  INTEGER :: ft
  INTEGER :: icode          ! Error return code (=0 is OK)
  INTEGER :: info           ! Dummy arg for GCOM calls
  INTEGER :: oasis_info     ! Return code from prism calls
  INTEGER :: infosst, infou, infov, infofreeze
  INTEGER :: infofreezen, infosnowthickn, infohicen
  INTEGER :: infomean
  INTEGER :: infoco2, infoco2fx

! gol124: auscom coupling
  REAL :: ltfs
  INTEGER :: infosal

! measure time spent in different stages of geto2a:
! PREP - before communication
! COMM - receiving data from oasis
! MEAN_POLAR_ROW - meaning
! POST - post processing received data
  IF (LTIMER) THEN
! DEPENDS ON: timer
      CALL TIMER('GETO2A_PREP',3)
  END IF

  ltfs = access_tfs

  ! Initialise potential incoming fields
  ocn_sst(:,:) = 0.0
  ocn_u(:,:) = 0.0
  ocn_v(:,:) = 0.0
  ocn_freeze(:,:) = 0.0
  ocn_freezen(:,:,:) = 0.0
  ocn_hice(:,:) = 0.0
  ocn_hicen(:,:,:) = 0.0
  ocn_snowthick(:,:) = 0.0
  ocn_snowthickn(:,:,:) = 0.0
  ocn_co2(:,:) = 0.0
  ocn_co2fx(:,:) = 0.0

  tstar_local(:,:) = 0.0
  tstar_ssi(:,:) = 0.0
  tstar_sice_local(:,:) = 0.0
  tstar_land_local(:,:) = 0.0

! gol124: auscom coupling
!  auscom_salinity = 0.0
  IF (.not.ocn_sss) THEN
      auscom_salinity(:,:) = (access_tfs-ZeroDegC)/-0.054
  END IF

  ! This saves us having to do an explicit mask because
  ! OASIS4 will just overlay the sea points
  IF (L_CTILE) THEN
    DO J = 1, oasis_jmt
      DO I = 1, oasis_imt
        ocn_sst(I,J)=TSTAR_SEA(I,J)
        ocn_sst_orig(I,J)=TSTAR_SEA(I,J)
      END DO
    END DO
  ELSE
    DO J = 1, oasis_jmt
      DO I = 1, oasis_imt
        ocn_sst(I,J)=TSTAR(I,J)
        ocn_sst_orig(I,J)=TSTAR(I,J)
      END DO
    END DO
  END IF


  DO J = 1, oasis_jmt
    DO I = 1, oasis_imt
      tstar_sice_local(I,J)=TSTAR_SICE(I,J)
      tstar_land_local(I,J)=TSTAR_LAND(I,J)
    END DO
  END DO

  DO K=1,nice
    DO J = 1, oasis_jmt
      DO I = 1, oasis_imt
        ocn_freezen(I,J,K)=ICE_FRACT_CAT(I,J,K)
        ocn_hicen(I,J,K)=ICE_THICK_CAT(I,J,K)
        ocn_snowthickn(I,J,K)=SNODEP_SEA_CAT(I,J,K)
      END DO
    END DO
  END DO

  IF (LTIMER) THEN
! DEPENDS ON: timer
      CALL TIMER('GETO2A_PREP',4)
  END IF

  IF (GET_STEP) THEN

    IF (LTIMER) THEN
! DEPENDS ON: timer
       CALL TIMER('GETO2A_COMM',3)
    END IF

    ft = fld_type_p

    IF (L_COUPLE_MASTER) THEN
      r_cols = glsize(1,ft)
      r_rows = glsize(2,ft)
    ELSE
      r_cols = oasis_imt
      r_rows = oasis_jmt
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_sst"
    var_ind = vind_ocn_sst
    ! DEPENDS ON: oasis3_get64
    CALL oasis3_get64(ocn_sst,oasis_imt,oasis_jmt,infosst,ft          &
       ,r_cols,r_rows,get_step                          &
       ,halo_type_no_halo,gc_all_proc_group,mype)
    IF (rdump_enable) THEN
        CALL write_rdata (dump_sst,ft,ocn_sst,oasis_imt,oasis_jmt,.false.)
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_freeze 1-5"
    infomean = 0
    DO K = 1, nice  ! Note: NICE must=5
      var_ind = vind_ocn_freezen(K)

      ! DEPENDS ON: oasis3_get64
      CALL oasis3_get64(ocn_freezen(1,1,k),oasis_imt,oasis_jmt       &
         ,infofreezen,ft                                  &
         ,r_cols,r_rows,get_step              &
         ,halo_type_no_halo,gc_all_proc_group,mype)
      IF (rdump_enable) THEN
          CALL write_rdata (dump_frzn(k), ft, ocn_freezen(1,1,k),    &
                          oasis_imt, oasis_jmt, .false.)
    END IF

      IF (K >  1.AND.infofreezen /= infomean) THEN
        WRITE(6,*) "INFORFREEZEN inconsistent"
        STOP
      END IF
      infomean=infofreezen
    END DO ! Over K

    write(6,*) "oasis3_geto2a: Get ocn snowthick 1-5"
    infomean = 0
    DO K = 1, nice  ! Note: NICE must=5
      var_ind = vind_ocn_snowthickn(K)

      ! DEPENDS ON: oasis3_get64
      CALL oasis3_get64(ocn_snowthickn(1,1,k),oasis_imt,oasis_jmt    &
         ,infosnowthickn,ft                               &
         ,r_cols,r_rows,get_step              &
         ,halo_type_no_halo,gc_all_proc_group,mype)
      IF (rdump_enable) THEN
          CALL write_rdata (dump_snwtn(k), ft, ocn_snowthickn(1,1,k),&
                            oasis_imt, oasis_jmt, .false.)
      END IF
      IF (K >  1.AND.infosnowthickn /= infomean) THEN
        WRITE(6,*) "INFOSNOWTHICKN inconsistent"
        STOP
      END IF
      infomean=infosnowthickn
    END DO

    write(6,*) "oasis3_geto2a: Get ocn_hice 1-5"
    infomean = 0
    DO K = 1, nice  ! Note: NICE must=5
      var_ind = vind_ocn_hicen(K)

      ! DEPENDS ON: oasis3_get64
      CALL oasis3_get64(ocn_hicen(1,1,k),oasis_imt,oasis_jmt         &
         ,infohicen,ft                                    &
         ,r_cols,r_rows,get_step              &
         ,halo_type_no_halo,gc_all_proc_group,mype)
      IF (rdump_enable) THEN
          CALL write_rdata (dump_hicn(k), ft, ocn_hicen(1,1,k),      &
                            oasis_imt, oasis_jmt, .false.)
      END IF
      IF (K >  1.AND.infohicen /= infomean) THEN
        WRITE(6,*) "INFOHICEN inconsistent"
        STOP
      END IF
      infomean=infohicen
    END DO

    ft = fld_type_u
    IF (L_COUPLE_MASTER) THEN
      r_cols = glsize(1,ft)
      r_rows = glsize(2,ft)
    ELSE
      r_cols = oasis_imt
      r_rows = oasis_jmt_u
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_u"
    var_ind = vind_ocn_u

    ! DEPENDS ON: oasis3_get64
    CALL oasis3_get64(ocn_u,oasis_imt,oasis_jmt_u,infou,ft            &
       ,r_cols,r_rows,get_step                          &
       ,halo_type_no_halo,gc_all_proc_group,mype)
    IF (rdump_enable) THEN
        CALL write_rdata (dump_suno, ft, ocn_u,              &
                          oasis_imt, oasis_jmt, .false.)
    END IF

    ft = fld_type_v
    IF (L_COUPLE_MASTER) THEN
      r_cols = glsize(1,ft)
      r_rows = glsize(2,ft)
    ELSE
      r_cols = oasis_imt
      r_rows = oasis_jmt_v
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_v"
    var_ind = vind_ocn_v

    ! DEPENDS ON: oasis3_get64
    CALL oasis3_get64(ocn_v,oasis_imt,oasis_jmt_v,infov,ft            &
       ,r_cols,r_rows,get_step                          &
       ,halo_type_no_halo,gc_all_proc_group,mype)
    IF (rdump_enable) THEN
        CALL write_rdata (dump_svno, ft, ocn_v,        &
                          oasis_imt, oasis_jmt, .true.)
    END IF

    ft = fld_type_p

    IF (L_COUPLE_MASTER) THEN
      r_cols = glsize(1,ft)
      r_rows = glsize(2,ft)
    ELSE
      r_cols = oasis_imt
      r_rows = oasis_jmt
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_co2"
    var_ind = vind_ocn_co2
    ! DEPENDS ON: oasis3_get64
    CALL oasis3_get64(ocn_co2,oasis_imt,oasis_jmt,infoco2,ft          &
       ,r_cols,r_rows,get_step                          &
       ,halo_type_no_halo,gc_all_proc_group,mype)
    IF (rdump_enable) THEN
        CALL write_rdata (dump_ocn_co2,ft,ocn_co2,oasis_imt,oasis_jmt,.false.)
    END IF
    
    !the oasis3_get ocn_co2 interface is here, but received ocn_co2 is not used in atmosphere


    ft = fld_type_p

    IF (L_COUPLE_MASTER) THEN
      r_cols = glsize(1,ft)
      r_rows = glsize(2,ft)
    ELSE
      r_cols = oasis_imt
      r_rows = oasis_jmt
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_co2fx"
    var_ind = vind_ocn_co2fx
    ! DEPENDS ON: oasis3_get64
    CALL oasis3_get64(ocn_co2fx,oasis_imt,oasis_jmt,infoco2fx,ft          &
       ,r_cols,r_rows,get_step                          &
       ,halo_type_no_halo,gc_all_proc_group,mype)
    IF (rdump_enable) THEN
        CALL write_rdata (dump_ocn_co2fx,ft,ocn_co2fx,oasis_imt,oasis_jmt,.false.)
    END IF
    !ocn_co2flux is assigned to atmosphere
    co2flux=ocn_co2fx

! gol124: auscom coupling
    IF (ocn_sss) THEN
        ft = fld_type_p

        IF (L_COUPLE_MASTER) THEN
            r_cols = glsize(1,ft)
            r_rows = glsize(2,ft)
        ELSE
            r_cols = oasis_imt
            r_rows = oasis_jmt
        END IF

        write(6,*) "oasis3_geto2a: Get ocn_sal"
        var_ind = vind_sal
! DEPENDS ON: oasis3_get64
        CALL oasis3_get64(auscom_salinity,oasis_imt,oasis_jmt,infosal,ft  &
                         ,r_cols,r_rows,get_step                          &
                         ,halo_type_no_halo,gc_all_proc_group,mype)
        IF (rdump_enable) THEN
            CALL write_rdata (dump_sss, ft, auscom_salinity,         &
                              oasis_imt, oasis_jmt, .false.)
        END IF
    ELSE
        infosal = infosst
    END IF

    IF (LTIMER) THEN
! DEPENDS ON: timer
       CALL TIMER('GETO2A_COMM',4)
    END IF

    IF (rdump_enable) THEN
        ! gol124: inc record counter, to be called after last write only!
        CALL inc_rrec
    END IF

    ! If this is a genuine get step.
    ! Should really be able to use the PRISM info
    ! flags to control this entirely, but it all adds to
    ! the run time.

    IF (L_COUPLE_MASTER) THEN
      ! Ensure all PEs know the return codes under all
      ! conditions - PE 0 should definitely have them
      CALL GC_IBCAST(805,1,0,nproc,info,infosst)
      CALL GC_IBCAST(810,1,0,nproc,info,infofreezen)
      CALL GC_IBCAST(825,1,0,nproc,info,infosnowthickn)
      CALL GC_IBCAST(830,1,0,nproc,info,infohicen)
      CALL GC_IBCAST(815,1,0,nproc,info,infou)
      CALL GC_IBCAST(820,1,0,nproc,info,infov)
! gol124: auscom coupling
      IF (ocn_sss) THEN
          CALL GC_IBCAST(835,1,0,nproc,info,infosal)
      ELSE
          infosal = infosst
      END IF
    END IF

    ! Check for successful receiving of fields from all possible
    ! sources.
    IF ((infosst == infofreezen).AND.                                 &
       (infosst == infohicen).AND.                                   &
       (infosst == infosnowthickn).AND.                              &
       (infosst == infosal).AND.                                     &
       (infosst == PRISM_Recvd.OR.                                   &
       infosst == PRISM_FromRest.OR.                                &
       infosst == PRISM_RecvOut.OR.                                 &
       infosst == PRISM_FromRestOut)) THEN

      ! If this PE deals with the North polar row we must
      ! set a uniform temperature in all grid points along that
      ! row to prevent crashes in the atmos dynamics.
      ! That's all jolly well in a 1xn decomposition, but a
      ! nx1 or nxm composition means we have to do some gathering
      ! and scattering.
      ! There shouldnt be any implications for conservation.
      ! Note: MEAN_POLAR_ROW must be called for all T grid point
      ! fields which we receive.

      ! This is a special adjustment to get sensible
      ! numbers in the polar row. One might regard this as a fix
      ! to make up for shortcomings in oasis3/4 remapping where
      ! it fails to deal with the N pole as a singularity.
      IF (at_extremity(Pnorth)) THEN
        DO I = 1, oasis_imt
          ocn_sst(i,oasis_jmt)=ocn_sst(i,oasis_jmt-1)
          ocn_freezen(i,oasis_jmt,:)=ocn_freezen(i,oasis_jmt-1,:)
          ocn_hicen(i,oasis_jmt,:)=ocn_hicen(i,oasis_jmt-1,:)
          ocn_snowthickn(i,oasis_jmt,:)=                           &
             ocn_snowthickn(i,oasis_jmt-1,:)
        END DO
! gol124: auscom coupling
        IF (ocn_sss) THEN
            DO I = 1, oasis_imt
               auscom_salinity(i,oasis_jmt)=auscom_salinity(i,oasis_jmt-1)
            END DO
        END IF
      END IF

! measure time spent in mean_polar_row as whole group
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('MEAN_POLAR_ROW',3)
      END IF

      ! DEPENDS ON: mean_polar_row
      CALL MEAN_POLAR_ROW(OCN_SST)
! gol124: auscom coupling
      IF (ocn_sss) THEN
      ! DEPENDS ON: mean_polar_row
          CALL MEAN_POLAR_ROW(auscom_salinity)
      END IF

      DO K=1,nice
        ! DEPENDS ON: mean_polar_row
        CALL MEAN_POLAR_ROW(OCN_FREEZEN(1,1,K))
        ! DEPENDS ON: mean_polar_row
        CALL MEAN_POLAR_ROW(OCN_HICEN(1,1,K))
        ! DEPENDS ON: mean_polar_row
        CALL MEAN_POLAR_ROW(OCN_SNOWTHICKN(1,1,K))
      END DO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('MEAN_POLAR_ROW',4)
      END IF

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('GETO2A_POST',3)
      END IF

      ! Perform necessary processing as described, for example,
      ! in 'Operations 5' on  HadGEM3 coupling diagram. (See HadGEM3 
      ! coupling documentation)

      DO K=1,nice
        DO J = 1, oasis_jmt
          DO I = 1, oasis_imt

            ! Ensure ice-depth and snow thickness are not negative (which 
            ! they can be if 2nd order conservative coupling is being used).
            ! No need to do this here for ice fraction because any
            ! negative values will get picked up by the minimum 
            ! ice fraction check below. 

            ocn_hicen(I,J,K)=MAX(0.0,ocn_hicen(I,J,K))  
            ocn_snowthickn(I,J,K)=MAX(0.0,ocn_snowthickn(I,J,K)) 

            ! Here we ensure that regardless of what has happened to the ice 
            ! fraction and GBM ice depth during the coupling,
            ! the minimum local ice depth is 1cm (as in CICE). This is done by
            ! squashing the ice to a smaller area whilst conserving volume (so
            ! the GBM ice depth is unchanged).
            !
            ! Note that freezen is unchanged unless the local ice depth (i.e.
            ! hicen/freezen is less than 1.0e-02, in which case freezen is
            ! decreased so that hicen/freezen=1.0e-02

            ocn_freezen(I,J,K)=MIN(ocn_freezen(I,J,K),            &
               (ocn_hicen(I,J,K)/1.0e-02))

            ! Also impose a minimum ice fraction of 2.0e-04

            IF (ocn_freezen(I,J,K).LT.2.0e-04) THEN
              ocn_freezen(I,J,K)=0.0
              ocn_hicen(I,J,K)=0.0
              ocn_snowthickn(I,J,K)=0.0
            END IF

            ocn_freeze(I,J)=ocn_freeze(I,J)+ocn_freezen(I,J,K)
            ocn_hicen(I,J,K)=ocn_hicen(I,J,K)+                    &
               ocn_snowthickn(I,J,K)*(kappai/kappas)
            ocn_hice(I,J)=ocn_hice(I,J)+ocn_hicen(I,J,K)
            ocn_snowthick(I,J)=ocn_snowthick(I,J)+                &
               ocn_snowthickn(I,J,K)

            IF (ocn_freezen(I,J,K) >  0.0) THEN
              ocn_hicen(I,J,K)=ocn_hicen(I,J,K)/                 &
                 ocn_freezen(I,J,K)
              ocn_snowthickn(I,J,K)=ocn_snowthickn(I,J,K)*       &
                 rhosnow/ocn_freezen(I,J,K)

            END IF
          END DO
        END DO
      END DO

      DO J = 1, oasis_jmt
        DO I = 1, oasis_imt
          IF (ocn_freeze(I,J) >  0.0) THEN
            ocn_freeze(I,J)=MIN(ocn_freeze(I,J),1.0)
            ocn_snowthick(I,J)=ocn_snowthick(I,J)*                &
               rhosnow/ocn_freeze(I,J)
            ocn_hice(I,J)=ocn_hice(I,J)/ocn_freeze(I,J)
          END IF
        END DO
      END DO

      ! The incoming SST is expected in K NOT degrees C
      ! Any necessary transforms will be done in the PSMILE
      ! controlled via the SMIOC, namcouple or in the sending
      ! model component.
      DO J = 1, oasis_jmt
        DO I = 1, oasis_imt

          ! Special form of masking for T points - if
          ! value returned from coupling is approx abs zero
          ! we assume this is an unchanged point masked during
          ! the coupling process and we reinstate
          ! our original value.
          IF (ocn_sst(I,J) <  1.0) THEN
            ocn_sst(I,J)=ocn_sst_orig(I,J)
! gol124: auscom coupling
! for salinity enforce the same masking as for sea surface temp
            auscom_salinity(i,j) = 0.0
          END IF
        END DO
      END DO

      ! Process surface temperature fields here (as used to be done
      ! in SWAPA20 pre UM version 7.0).

      DO J = 1, oasis_jmt
        DO I = 1, oasis_imt
          tstar_ssi(I,J)=ocn_sst(I,J)

          IF(fland_loc(I,J) <  1.0) THEN
            IF(LTLEADS.EQV..FALSE.) THEN
               IF (ocn_sss) THEN
                   LTFS = ZeroDegC - 0.054 * auscom_salinity(I,J)
               END IF
              IF(ocn_freeze(I,J) >  0.0) ocn_sst(I,J)=LTFS
            END IF
            tstar_sice_local(I,J)=AMIN1(tstar_sice_local(I,J),TM)
            tstar_ssi(I,J)=ocn_freeze(I,J)*tstar_sice_local(I,J)+  &
               (1.0-ocn_freeze(I,J))*ocn_sst(I,J)
          END IF
          tstar_local(I,J)=fland_loc(I,J)*tstar_land_local(I,J)+    &
             (1.0-fland_loc(I,J))*tstar_ssi(I,J)
        END DO
      END DO

      ! Having got our new values we need to move them
      ! to the appropriate arrays where they'll be employed next 
      ! time the appropriate calculation is called.
      DO J = 1, oasis_jmt
        DO I = 1, oasis_imt

          IF (L_CTILE) THEN
            TSTAR_SEA(I,J)=ocn_sst(I,J)
          END IF

          TSTAR(I,J)=tstar_local(I,J)

          TSTAR_SICE(I,J)=tstar_sice_local(I,J)
          ICE_FRACTION(I,J,1)=ocn_freeze(I,J)
          ICE_THICKNESS(I,J,1)=ocn_hice(I,J)
          SNODEP_SEA(I,J)=ocn_snowthick(I,J)
        END DO
      END DO

      DO K = 1,nice
        DO J = 1, oasis_jmt
          DO I = 1, oasis_imt
            ICE_FRACT_CAT(I,J,K)=ocn_freezen(I,J,K)
            ICE_THICK_CAT(I,J,K)=ocn_hicen(I,J,K)
            SNODEP_SEA_CAT(I,J,K)=ocn_snowthickn(I,J,K)
          END DO
        END DO
      END DO

      IF (LTIMER) THEN
! DEPENDS ON: timer
          CALL TIMER('GETO2A_POST',4)
      END IF

    END IF

    ! Update our surface currents (U)
    IF ( infou == PRISM_Recvd        .OR.                             &
       infou == PRISM_FromRest     .OR.                             &
       infou == PRISM_RecvOut      .OR.                             &
       infou == PRISM_FromRestOut ) THEN

      ! Incoming currents expected in m/s

      ! This is a work around to get sensible
      ! numbers in the polar row -  it is NOT necessarily
      ! sensible for all configurations (e.g. LAMs). 
      ! It is a fix to make up for shortcomings in 
      ! remapping 
      IF (at_extremity(Pnorth)) THEN
        DO I = 1, oasis_imt
          ocn_u(i,oasis_jmt_u)=0.0
        END DO
      END IF

      DO J = 1, oasis_jmt_u
        DO I = 1, oasis_imt
          U_SEA(I,J)=ocn_u(I,J)
        END DO
      END DO
    END IF

    ! Update our surface currents (V)
    IF ( infov == PRISM_Recvd        .OR.                             &
       infov == PRISM_FromRest     .OR.                             &
       infov == PRISM_RecvOut      .OR.                             &
       infov == PRISM_FromRestOut ) THEN

      ! Incoming currents expected in m/s
      DO J = 1, oasis_jmt_v
        DO I = 1, oasis_imt
          V_SEA(I,J)=ocn_v(I,J)
        END DO
      END DO
    END IF

  END IF ! If GET_STEP = true

END SUBROUTINE oasis3_geto2a
