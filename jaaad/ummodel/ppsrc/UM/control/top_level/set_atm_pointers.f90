
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SET_ATM_POINTERS ---------------------------------------
!LL
!LL  Set pointers for primary atmosphere fields
!LL Initialisation routine for CRAY YMP
!LL
!LL MC, CW      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.1    9/02/93  : added comdeck CHSUNITS to define NUNITS for
!LL                    comdeck CCONTROL.
!LL 3.1     11/03/93  Set JTRACER(1,1) to have a sensible address
!L                    even if there are no tracers to remove bounds
!L                    checking problems in later routines. R. Rawlins
!LL 3.2    27/03/93 Dynamic allocation of main data arrays. R. Rawlins
!LL 3.4    20/01/94 Changes to allow for non-consecutive tracers.
!LL                 M.Carter
!LL  3.3  22/11/93  Add aerosol ancillary fields.  R T H Barnes.
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                   Remove a_max_variables.
!LL                   PPINDEX now read from UI in INITCTL.
!LL  3.4  05/09/94  Add murk & user ancillary fields.  RTHBarnes.
!LL  3.4   18/05/94  J.Thomson add pointers for slab temperature
!LL                            and u and v ice velocities.
!LL  3.5   19/05/95  Some pointers for level dependent constants
!LL                  array removed. Sub_model change. D. Robinson
!LL  4.0   06/09/95  Set up pointers correctly for SLAB fields.
!LL                  D. Robinson
!LL
!LL  4.0 26/7/94 R.E. Essery Extra prognostic sea-ice temperature.
!LL  4.1 10/1/96      Extra prognostics frozen and unfrozen soil
!LL                   moisture fractions, and canopy conductance
!LL                   plus 2 extra vegetation fields J.Smith
!LL  4.1 30/04/96     Add pointers for 6 new variables and 6 new
!LL                   ancillary fields for Sulphur Cycle   M Woodage
!LL  4.3  18/3/97     Add pointers for HadCM2 sulphate loading patterns
!LL                                                   William Ingram
!LL  4.2  16/08/96    Added mpp PARVARS comdeck, and defined rowdepc
!LL                   (filterable wave numbers array) to be globally
!LL                   sized.                                P.Burton
!LL  4.3  06/03/97    Dimension multi_level land fields by LAND_FIELD
!LL                   for mpp jobs. D. Robinson.
!LL  4.4  05/09/97    Add pointer for net energy flux prognostic
!LL                   S.D.Mullerworth
!LL  4.4  04/08/97  Generalise JQCF pointer for mixed phase
!LL                 precipitation scheme.  RTHBarnes.
!LL  4.4  05/08/97    Add pointer for convective cloud amount on
!LL                   model levels (3D CCA) if L_3D_CCA. J.M.Gregory
!LL  4.4  10/09/97    Added pointers for snow grain size and snow soot
!LL                   content used in prognostic snow albedo scheme
!LL                   R. Essery
!LL  4.4  16/09/97    Add call to NSTYPES and pointers for new
!LL                   vegetation and land surface prognostics. R.A.Betts
!LL  4.5   1/07/98    Add pointers for ocean CO2 flux and surface
!LL                     CO2 emissions. C.D.Jones
!LL  4.5  04/03/98   Add pointers for NH3 prognostic and NH3 surface
!LL                    emiss for S Cycle                 M Woodage
!LL                  Also add pointers for 3 soot prognostic variables
!LL                    and 2 soot emiss                    M Woodage
!LL  4.5  15/07/98  Add pointers for new 3D CO2 array. C.D.Jones
!    4.5  22/10/98    Add pointers for extra multi-layer user
!                     ancillary fields
!                     Author D.M. Goddard
!    4.5  29/04/98    Pointer to total soil moisture content to point
!                     to non prognostic space in MOSES.
!                     Author D.M. Goddard
!LL  4.5  19/01/98    Replace JVEG_FLDS and JSOIL_FLDS with
!LL                   individual pointers. D. Robinson
!LL  5.0  13/07/99    Extensive changes for C-P C grid upgrade
!LL                   D Goddard
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!LL  5.1  28/02/00    Add pointers for ND LBCs.
!LL                   Old JRIM pointers are kept in until all
!LL                   references in code are deleted at 5.2
!LL                                                         P.Burton
!  5.1  28/04/00   Extra pressure level above top theta level
!                  p_top_theta_level removed                A.Malcolm
!LL  5.2  13/09/00    Remove JSOOT, JSO4, JH2SO4 pointers and
!LL                   change sizes for tracers. P.Selwood
!LL
!LL  5.1  29/02/00    Add JNET_MFLUX R A Stratton.
!LL  5.2  18/10/00    Removed JUICE and JVICE.            K.D.Williams
!    5.2  31/08/00    Add JOROG_LBC pointer
!                     Remove levels from LBC pointers
!                     Change section/item for LBC pointers
!                     Remove JRIM and JRIM_TENDENCY
!                                                   P.Burton
!    5.2  25/09/00    Clear out dead legacy code (3D RHcrit not needed
!                     in dump).                           A.C.Bushell
!
!LL  5.2  15/11/00  Additional initialisation of pointers for
!LL                 MOSES 2 and dynamic vegetation         M. Best
!LL  5.2  09/03/01  Introduce extra pointers in level dependent
!                   constants for height definitions. R Rawlins
!LL  5.3  04/10/01  Added new Slab model prognostics and re-order above
!LL                 where Sect_No is set to 31.            K.Williams
!    5.3  20/07/01  Replace hardwires with access of SI index when
!                   calculating pointers for items 406 to 412
!                   S.D.Mullerworth
!LL  5.3  19/06/01   Add JTPPSOZONE Dave Tan
!    5.4  02/09/02  New pointer JSNODEP_SEA added        K.Williams
!    5.4  28/08/02  Set up pointers for canopy snow and snow
!                   beneath canopy.  R. Essery
!    5.4  27/07/02  Increase cloud fraction halo sizes for PC2
!                                                          D. Wilson
!    5.5  05/11/02  Large-scale hydrology scheme: Add pointers for
!                   gamma function, water table depth,
!                   surface saturated fraction and wetland fraction.
!                                                        N. Gedney
!    5.5  05/02/03  Add pointers for biomass aerosol scheme   P Davison
!    5.5  13/02/03  New pointers for multicatagory sea ice. J.Ridley
!
!    5.5  03/02/03  Add pointers for 3 moisture arrays
!                   (qcf2,qrain,qgraup) and lbc,lbc_tend     R.M.Forbes
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
!
!    5.5  26/02/03  Add pointers for river routing. P.Selwood.
!  6.0  14/08/03  NEC SX-6 Fix - avoid undefined value?
!                 R Barnes & J-C Rioual.
!    6.0  30/07/03  Add pointers for cloud fraction lbcs. Damian Wilson
!  6.0    12/08/03 removes JRIV_GRIDBOX. C.Bunton.
!    6.2  11/11/05  Set pointers for section 34 UKCA tracers. R Barnes.
!    6.1  28/06/04  Set pointers for section 33 free tracers. R Barnes.
!  6.1    20/08/03  Add pointers for STOCHEM fields.  C. Johnson
!    6.1  07/04/04  Add pointers for DMS concentration.      A. Jones
!    6.2  13/07/05  Correct free tracer LBC pointers. R Barnes.
!    6.2  19/08/05  Fix out of bounds reference. P.Selwood.
!    6.2  15/08/05  Free format fixes. P.Selwood
! 6.2  02/03/06 Added pointer direct PAR flux. M.G. Sanderson
!    6.2  01/10/04  Add pointers for murk aerosol lbc.    R.M.Forbes
!    6.2  01/03/06  Add pointers for RothC pools and fluxes. C.D. Jones
!   6.2   21/2/06  Re-route outflow from inland basins to soil moisture
!                  P. Falloon

! 6.2  25/11/05 Functionality for improved time stepping, radiative
!               forcing and radiance code added for versions 3C
!               and 3Z of radiation code             (J.-C. Thelen)
!    6.2  24/02/06  Add pointers for ocean DMS flux             J.Gunson
!LL Programming Standard: Unified Model DP NO. 3, Version 3
!LL
!LL  Logical task: P0
!LL
!LL  System Components: C21 (Atmosphere part)
!LL
!LL  Purpose:   Sets integer pointers to atmospheric
!LL             variables from STASHIN addresses.
!LL
!LL  External documentation: UMDP NO. C4 Version NO. 4
!LL
!LLEND-------------------------------------------------------------

      SUBROUTINE SET_ATM_POINTERS(                                      &
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
! ARGPTRA start
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add murk and user ancillary pointers. RTHBarnes
!  4.1  04/12/95  Add pointers JSTHU and JSTHF. J.Smith
!  4.1  26/04/96  Add pointers for Sulphur Cycle variables (12)  MJW
!  4.3   18/3/97  And for HadCM2 sulphate loading patterns.  Will Ingram
!  4.4   05/8/97  And for Conv. cloud amt on model levs. Julie Gregory
!  4.5  04/03/98   Remove pointer SO2_HILEM (add to CARGPT_ATMOS)
!                  Add 1 pointers for NH3 in S Cycle
!                  Add 3 pointers for Soot              M. Woodage
!  4.5  08/05/98   Add 16 new pointers for User Anc.    D. Goddard
!  4.5  15/07/98   Add pointers for new 3D CO2 array.   C.D.Jones
!  4.5  17/08/98   Remove JSOIL_FLDS and JVEG_FLDS      D. Robinson
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  5.1  28/02/00   Added ND LBC variables            P.Burton
!  5.2  31/08/00   Added orography LBC and changed order of LBC pointers
!                                                               P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2    13/09/00 Remove JSOOT, JSO4, JH2SO4 pointers. P.Selwood.
!  5.3  19/06/01   Add pointer (JTPPSOZONE) for
!                  tropopause-based ozone        D. Tan
!  5.5  01/11/02   Add pointers JQCF2, JQRAIN, and JQGRAUP
!                  and associated LBCs for microphysics    R.M.Forbes
!  6.1  06/11/03   Add pointers for STOCHEM fields
!  6.0  30/07/03   Add pointers for cloud fractions and cloud
!                  fraction lbcs. Damian Wilson
! 6.1  13/08/04  Reduce number of continuation lines - prevent coupled
!                models hitting 99 line limit. R. Hill
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  11/05/05   Put STOCHEM pointers in correct place.
!                      M.G. Sanderson
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add JMURK_LBC pointers.    R.M.Forbes
      ! Pointers for ATMOSPHERE model variables. Configuration dependent.

      ! Addresses in D1 array of primary variables held in primary and
      ! secondary space

      !  Array  variables (dimensions are resolution dependent)

      ! Data variables stored in primary space.
     &  JU, JV, JW, JRHO, JTHETA, JQ, JQCL, JQCF,                       &
     &  JEXNER_RHO_LEVELS, JU_ADV, JV_ADV, JW_ADV,                      &
      ! Data variables stored in secondary space.
     &  JP, JP_THETA_LEVELS,  JEXNER_THETA_LEVELS,                      &
      ! Cloud Fields
     &  JCCA, JCF_AREA, JCF_BULK, JCF_LIQUID, JCF_FROZEN,               &
      ! Soil Ancillary Fields
     &  J_DEEP_SOIL_TEMP,  JSMCL, JSTHU, JSTHF,                         &
      ! Radiation increments
     &  JSW_INCS, JLW_INCS,                                             &
      ! Ozone
     &  JOZONE,                                                         &
      ! Tracers and Aerosols
     &  JTRACER, JMURK, JMURK_SOURCE,                                   &
     &  JSO2, JDMS, JSO4_AITKEN, JSO4_ACCU, JSO4_DISS, JH2O2,           &
     &  JNH3, JSOOT_NEW, JSOOT_AGD, JSOOT_CLD, JSO2_NATEM,              &
     &  JOH, JHO2, JH2O2_LIMIT,JO3_CHEM, JHadCM2_SO4, JCO2,             &
      ! User Ancillary fields
     &  JUSER_MULT1, JUSER_MULT2, JUSER_MULT3, JUSER_MULT4,             &
     &  JUSER_MULT5, JUSER_MULT6, JUSER_MULT7, JUSER_MULT8,             &
     &  JUSER_MULT9, JUSER_MULT10, JUSER_MULT11, JUSER_MULT12,          &
     &  JUSER_MULT13, JUSER_MULT14, JUSER_MULT15, JUSER_MULT16,         &
     &  JUSER_MULT17, JUSER_MULT18, JUSER_MULT19, JUSER_MULT20,         &
      ! Lateral Boundary Conditions and tendencies
     &  JOROG_LBC,JU_LBC,JU_LBC_TEND,JV_LBC,JV_LBC_TEND,JW_LBC,         &
     &  JW_LBC_TEND,JRHO_LBC,JRHO_LBC_TEND,JTHETA_LBC,JTHETA_LBC_TEND,  &
     &  JQ_LBC,JQ_LBC_TEND,JQCL_LBC, JQCL_LBC_TEND,JQCF_LBC,            &
     &  JQCF_LBC_TEND,JEXNER_LBC,JEXNER_LBC_TEND,JU_ADV_LBC,            &
     &  JU_ADV_LBC_TEND,JV_ADV_LBC,JV_ADV_LBC_TEND,JW_ADV_LBC,          &
     &  JW_ADV_LBC_TEND,JTRACER_LBC,JTRACER_LBC_TEND,                   &
! Tropopause-based Ozone
     &  JTPPSOZONE,                                                     &
      ! Biomass aerosol
     &  JBMASS_NEW, JBMASS_AGD, JBMASS_CLD,                             &
      ! Additional microphysics fields and lbcs
     &  JQCF2,JQRAIN,JQGRAUP,JQCF2_LBC,JQCF2_LBC_TEND,JQRAIN_LBC,       &
     &  JQRAIN_LBC_TEND,JQGRAUP_LBC,JQGRAUP_LBC_TEND,JDUST_DIV1,        &
     &  JDUST_DIV2,JDUST_DIV3,                                          &
     &  JDUST_DIV4,JDUST_DIV5,JDUST_DIV6,                               &
     &  JCF_BULK_LBC,JCF_BULK_LBC_TEND,JCF_LIQUID_LBC,                  &
     &  JCF_LIQUID_LBC_TEND,JCF_FROZEN_LBC,JCF_FROZEN_LBC_TEND,         &
! Fields from STOCHEM
     &  JCH4_STOCH,JO3_STOCH,                                           &
! Pointer for direct PAR flux
     &  JDIRPAR,                                                        &
     &  JTR_UKCA, JMURK_LBC, JMURK_LBC_TEND,                            &
! Pointers for UKCA oxidant fields
     &  JOH_UKCA, JHO2_UKCA, JH2O2_UKCA, JO3_UKCA,                      &
! Convective Cloud Fields
     &  JLCBASE, JCCW_RAD,                                              &
! Ozone tracer and cariolle parameters
     &  JOZONE_TRACER,JO3_PROD_LOSS,JO3_P_L_VMR,JO3_VMR,JO3_P_L_TEMP,   &
     &  JO3_TEMP,JO3_P_L_COLO3,JO3_COLO3,                               &
! Pointers for Aerosol climatologies
     &  JARCLBIOG_BG, JARCLBIOM_FR, JARCLBIOM_AG, JARCLBIOM_IC,         &
     &  JARCLBLCK_FR, JARCLBLCK_AG, JARCLSSLT_FI, JARCLSSLT_JT,         &
     &  JARCLSULP_AC, JARCLSULP_AK, JARCLSULP_DI, JARCLDUST_B1,         &
     &  JARCLDUST_B2, JARCLDUST_B3, JARCLDUST_B4, JARCLDUST_B5,         & 
     &  JARCLDUST_B6, JARCLOCFF_FR, JARCLOCFF_AG, JARCLOCFF_IC,         &
     &  JARCLDLTA_DL,                                                   &
! Fossil-fuel organic carbon aerosol
     &  JOCFF_NEW, JOCFF_AGD, JOCFF_CLD,                                &
      ! CABLE
     &  JTSOIL_TILE, JSMCL_TILE, JSTHF_TILE, JSNOW_DEPTH3L,             &
     &  JSNOW_MASS3L, JSNOW_TMP3L, JSNOW_RHO3L, JSNOW_RHO1L, JSNOW_AGE, & 
     &  JSNOW_FLG3L,                                                    &
! Lestevens Sept 2012 - adding new progs
     &  JCPOOL_TILE,JNPOOL_TILE,JPPOOL_TILE,JSOIL_ORDER,                &
     &  JNIDEP,JNIFIX,JPWEA,JPDUST,JGLAI,JPHENPHASE,                    &
     &  JCABLE_LAI,                                                     &
! ARGPTRA end
     &                  ICODE,CMESSAGE)

      USE CSENARIO_MOD
      IMPLICIT NONE

!L
!*L Arguments
!L
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
      INTEGER                                                           &
     &    ICODE                  ! OUT: Error return code
!
      CHARACTER*80                                                      &
     &    CMESSAGE               ! OUT: Error return message
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

!*L------------------ COMDECK LOOKADD ----------------------------------
!LL
!LL Purpose : Contains information about the format
!LL           of the PP header
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
!LL                   BRSVD2 to BHULEV and definitions for BRLEV and
!LL                   BHRLEV. Corresponding changes made to STWORK1A
!LL                   and PPHEAD1A. (Andrew Brady)
!LL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
!LL  5.1  17/04/00    Fixed/Free format. P.Selwood.
!LL  5.2  25/09/00    Add LBCC_xxxx variables for the compressed
!LL                   LBC LOOKUP array                  P.Burton
!LL
!LL Programming standard :
!LL
!LL Logical components covered : F092
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LLEND -----------------------------------------------------------------
!

! Validity time
      Integer, Parameter :: LBYR   =1   ! Year
      Integer, Parameter :: LBMON  =2   ! Month
      Integer, Parameter :: LBDAT  =3   ! Day of month
      Integer, Parameter :: LBHR   =4   ! Hour
      Integer, Parameter :: LBMIN  =5   ! Minute
      Integer, Parameter :: LBDAY  =6   ! Day number

! Data time
      Integer, Parameter :: LBYRD  =7   ! Year
      Integer, Parameter :: LBMOND =8   ! Month
      Integer, Parameter :: LBDATD =9   ! Day of month
      Integer, Parameter :: LBHRD  =10  ! Hour
      Integer, Parameter :: LBMIND =11  ! Minute
      Integer, Parameter :: LBDAYD =12  ! Day number

      Integer, Parameter :: LBTIM  =13  ! Time indicator
      Integer, Parameter :: LBFT   =14  ! Forcast period (hours)
      Integer, Parameter :: LBLREC =15  ! Length of data record
      Integer, Parameter :: LBCODE =16  ! Grid type code
      Integer, Parameter :: LBHEM  =17  ! Hemisphere indicator
      Integer, Parameter :: LBROW  =18  ! Number of rows in grid
      Integer, Parameter :: LBNPT  =19  ! Number of points per row
      Integer, Parameter :: LBEXT  =20  ! Length of extra data
      Integer, Parameter :: LBPACK =21  ! Packing method indicator
      Integer, Parameter :: LBREL  =22  ! Header release number
      Integer, Parameter :: LBFC   =23  ! Field code
      Integer, Parameter :: LBCFC  =24  ! Second field code
      Integer, Parameter :: LBPROC =25  ! Processing code
      Integer, Parameter :: LBVC   =26  ! Vertical coordinate type
      Integer, Parameter :: LBRVC  =27  ! Coordinate type for reference
                                        ! level

      Integer, Parameter :: LBEXP  =28  ! Experiment number
      Integer, Parameter :: LBEGIN =29  ! Start record
      Integer, Parameter :: LBNREC =30  ! No of records-Direct access
                                        ! only
      Integer, Parameter :: LBPROJ =31  ! Met-O-8 projection number
      Integer, Parameter :: LBTYP  =32  ! Met-O-8 field type
      Integer, Parameter :: LBLEV  =33  ! Met-O-8 level code
      Integer, Parameter :: LBRSVD1=34  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD2=35  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD3=36  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD4=37  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBSRCE =38  ! =1111 to indicate following
                                        ! apply to UM
      Integer, Parameter :: DATA_TYPE =39  ! Indicator for real/int or
                                           ! timeseries
      Integer, Parameter :: NADDR  =40  ! Start address in DATA_REAL or
                                        ! DATA_INT
      Integer, Parameter :: LBUSER3=41  ! Free for user-defined function
      Integer, Parameter :: ITEM_CODE =42  !Stash code
      Integer, Parameter :: LBPLEV =43  ! Pseudo-level indicator (if
                                        ! defined)
      Integer, Parameter :: LBUSER6=44  ! Free for user-defined function
      Integer, Parameter :: MODEL_CODE =45 ! internal model identifier

      Integer, Parameter :: BULEV  =46  ! Upper level boundary
      Integer, Parameter :: BHULEV =47  ! Upper level boundary
      Integer, Parameter :: BRSVD3 =48  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BRSVD4 =49  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BDATUM =50  ! Datum value
      Integer, Parameter :: BACC   =51  ! (Packed fields) Packing
                                        ! accuracy
      Integer, Parameter :: BLEV   =52  ! Level
      Integer, Parameter :: BRLEV  =53  ! Lower level boundary
      Integer, Parameter :: BHLEV  =54  ! (Hybrid levels) A-level of
                                        ! value
      Integer, Parameter :: BHRLEV =55  ! Lower level boundary
      Integer, Parameter :: BPLAT  =56  ! Real latitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BPLON  =57  ! Real longitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BGOR   =58  ! Grid orientation
      Integer, Parameter :: BZY    =59  ! Zeroth latitude
      Integer, Parameter :: BDY    =60  ! Latitude interval
      Integer, Parameter :: BZX    =61  ! Zeroth longitude
      Integer, Parameter :: BDX    =62  ! Longitude interval
      Integer, Parameter :: BMDI   =63  ! Missing data indicator
      Integer, Parameter :: BMKS   =64  ! M,K,S scaling factor

      Integer, Parameter :: LBCC_LBYR   = 1  ! Year
      Integer, Parameter :: LBCC_LBMON  = 2  ! Month
      Integer, Parameter :: LBCC_LBDAT  = 3  ! Day of the month
      Integer, Parameter :: LBCC_LBHR   = 4  ! Hour
      Integer, Parameter :: LBCC_LBMIN  = 5  ! Minute
      Integer, Parameter :: LBCC_LBDAY  = 6  ! Day number
      Integer, Parameter :: LBCC_LBEGIN = 7  ! Start record
      Integer, Parameter :: LBCC_NADDR  = 8  ! Start address of DATA
! Mapping of MPP_LOOKUP; analogous to mapping in PP header

      Integer, Parameter :: P_NADDR=1    ! Address on local PE
      Integer, Parameter :: P_LBLREC=2   ! Local length of record

!*----------------------------------------------------------------------
! NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
! ITEM_CODE is the location in PP header for a code defined as
!           (section number)*1000+item number
! DATA_TYPE is the location in the PP header defining data as REAL or
!           INTEGER.
! LBNPT is the location defining the number of points per row
!
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

! local variables

      INTEGER                                                           &
     &        IVAR,                                                     &
                                 ! Loop counts
     &        JVAR,                                                     &
                                 ! Loop counts
     &        IFLD,                                                     &
     &        LEV                                                       &
     &       ,im_ident                                                  &
                            !  Internal Model Identifier
     &       ,im_index                                                  &
                            !  Internal Model Index in Stash arrays
     &       ,Sect_No       !  Stash section number

! Tracer fluxes stash section number - kdcorbin, 06/10
     INTEGER :: tflux_sect

!     Set to atmosphere internal model
      im_ident  = atmos_im
      im_index  = internal_model_index(im_ident)
      Sect_No   = 0


! Set pointers for atmospheric primary variables from STASH :-
      JU(1)                = SI(  2,Sect_No,im_index)
      JV(1)                = SI(  3,Sect_No,im_index)
      JTHETA(1)            = SI(  4,Sect_No,im_index)
      JQ(1)                = SI( 10,Sect_No,im_index)
      JQCF(1)              = SI( 12,Sect_No,im_index)
      JTSTAR               = SI( 24,Sect_No,im_index)
      JLAND                = SI( 30,Sect_No,im_index)
      JOROG                = SI( 33,Sect_No,im_index)
      JW(0)                = SI(150,Sect_No,im_index)
      JRHO(1)              = SI(253,Sect_No,im_index)
      JQCL(1)              = SI(254,Sect_No,im_index)
      JEXNER_RHO_LEVELS(1) = SI(255,Sect_No,im_index)
      JQCF2(1)             = SI(271,Sect_No,im_index)
      JQRAIN(1)            = SI(272,Sect_No,im_index)
      JQGRAUP(1)           = SI(273,Sect_No,im_index)

!LL  5.3     06/01  Introduce extra pointers for coastal tiling.
!                                                     Nic Gedney.
      JFRAC_LAND     = SI(505,Sect_No,im_index)
      JTSTAR_LAND    = SI(506,Sect_No,im_index)
      JTSTAR_SEA     = SI(507,Sect_No,im_index)
      JTSTAR_SICE    = SI(508,Sect_No,im_index)
! Set pointers for seaice and land albedos
      JSICE_ALB      = SI(509,Sect_No,im_index)
      JLAND_ALB      = SI(510,Sect_No,im_index)
!
!LL  5.3     06/01  Introduce extra pointers for large-scale hydrology.
!                                                     Nic Gedney.
      JTI_MEAN       = SI(274,Sect_No,im_index)
      JTI_SIG        = SI(275,Sect_No,im_index)
      JFEXP          = SI(276,Sect_No,im_index)
      JGAMMA_INT     = SI(277,Sect_No,im_index)
      JWATER_TABLE   = SI(278,Sect_No,im_index)
      JFSFC_SAT      = SI(279,Sect_No,im_index)
      JF_WETLAND     = SI(280,Sect_No,im_index)
      JSTHZW         = SI(281,Sect_No,im_index)
      JA_FSAT        = SI(282,Sect_No,im_index)
      JC_FSAT        = SI(283,Sect_No,im_index)
      JA_FWET        = SI(284,Sect_No,im_index)
      JC_FWET        = SI(285,Sect_No,im_index)
      DO LEV=1,MODEL_LEVELS
        JW(LEV)=JW(LEV-1)+theta_off_size
      END DO
      DO LEV=2,MODEL_LEVELS
        JU(LEV)    =JU(LEV-1)     + u_off_size
        JV(LEV)    =JV(LEV-1)     + v_off_size
        JRHO(LEV)  =JRHO(LEV-1)   + theta_off_size
        JTHETA(LEV)=JTHETA(LEV-1) + theta_off_size
      END DO
      DO LEV=2,MODEL_LEVELS+1
        JEXNER_RHO_LEVELS(LEV)=JEXNER_RHO_LEVELS(LEV-1)+theta_off_size
      END DO
      DO LEV=2,WET_LEVELS
        JQ(LEV)    =JQ(LEV-1)     + theta_halo_size
        JQCL(LEV)  =JQCL(LEV-1)   + theta_halo_size
        JQCF(LEV)  =JQCF(LEV-1)   + theta_halo_size
        JQCF2(LEV)   =JQCF2(LEV-1)   + theta_halo_size
        JQRAIN(LEV)  =JQRAIN(LEV-1)  + theta_halo_size
        JQGRAUP(LEV) =JQGRAUP(LEV-1) + theta_halo_size
      END DO

      ! Pointers required to save coupling fields as
      ! prognostics when employing OASIS as the coupler.
      ! In non-coupled models these pointers will
      ! simply end up with a value of 1.
      JC_SOLAR = SI(171,Sect_No,im_index)
      JC_BLUE =  SI(172,Sect_No,im_index)
      JC_DOWN =  SI(173,Sect_No,im_index)
      JC_LONGWAVE =  SI(174,Sect_No,im_index)
      JC_TAUX =  SI(176,Sect_No,im_index)
      JC_TAUY =  SI(177,Sect_No,im_index)
      JC_WINDMIX =  SI(178,Sect_No,im_index)
      JC_SENSIBLE =  SI(179,Sect_No,im_index)
      JC_SUBLIM =  SI(180,Sect_No,im_index)
      JC_EVAP =  SI(181,Sect_No,im_index)
      JC_BOTMELTN =  SI(184,Sect_No,im_index)
      JC_TOPMELTN =  SI(185,Sect_No,im_index)
      JC_LSRAIN =  SI(186,Sect_No,im_index)
      JC_LSSNOW =  SI(187,Sect_No,im_index)
      JC_CVRAIN =  SI(188,Sect_No,im_index)
      JC_CVSNOW =  SI(189,Sect_No,im_index)
      JC_RIVEROUT =  SI(192,Sect_No,im_index)

! gol124: auscom coupling
! (the same as JPSTAR)
      JC_PRESS = SI(409,Sect_No,im_index)

! Set pointers for optional atmospheric primary variables.
!  To be enabled at a later version than 5.0
!     IF(bit compatability required)
        JZH           = SI( 25,Sect_No,im_index)
        JU_ADV(1)     = SI(256,Sect_No,im_index)
        JV_ADV(1)     = SI(257,Sect_No,im_index)
        JW_ADV(0)     = SI(258,Sect_No,im_index)
        JNTML         = SI(259,Sect_No,im_index)
        JNBDSC        = SI(260,Sect_No,im_index)
        JNTDSC        = SI(261,Sect_No,im_index)
        JCUMULUS      = SI(262,Sect_No,im_index)
        JT1_SD        = SI(263,Sect_No,im_index)
        JQ1_SD        = SI(264,Sect_No,im_index)
        JCF_AREA(1)   = SI(265,Sect_No,im_index)
        JCF_BULK(1)   = SI(266,Sect_No,im_index)
        JCF_LIQUID(1) = SI(267,Sect_No,im_index)
        JCF_FROZEN(1) = SI(268,Sect_No,im_index)

        IF (L_3D_CCA) THEN
          JCCA(1) = SI(211,Sect_No,im_index)
          DO LEV=2,N_CCA_LEV
            JCCA(LEV)=JCCA(LEV-1)+THETA_FIELD_SIZE
          ENDDO
        ELSE
          JCCA(1)      = SI( 13,Sect_No,im_index)
        ENDIF

        JCCB           = SI( 14,Sect_No,im_index)
        JCCT           = SI( 15,Sect_No,im_index)
        JCCLWP         = SI( 16,Sect_No,im_index)
        JLCBASE        = SI( 21,Sect_No,im_index)
        JCANOPY_WATER  = SI( 22,Sect_No,im_index)
        JCCW_RAD(1)    = SI(212,Sect_No,im_index)

        DO LEV=2, WET_LEVELS
            JCCW_RAD(LEV)=JCCW_RAD(LEV-1)+THETA_FIELD_SIZE
        ENDDO
!     Endif !IF(bit compatability required)

      DO LEV=1,MODEL_LEVELS
        JW_ADV(LEV)=JW_ADV(LEV-1)+theta_halo_size
      END DO
      DO LEV=2,MODEL_LEVELS
        JU_ADV(LEV)=JU_ADV(LEV-1)+u_halo_size
        JV_ADV(LEV)=JV_ADV(LEV-1)+v_halo_size
      END DO
      DO LEV=2,WET_LEVELS
        JCF_AREA(LEV)  =JCF_AREA(LEV-1)+THETA_FIELD_SIZE
        JCF_BULK(LEV)  =JCF_BULK(LEV-1)+THETA_halo_SIZE
        JCF_LIQUID(LEV)=JCF_LIQUID(LEV-1)+THETA_halo_SIZE
        JCF_FROZEN(LEV)=JCF_FROZEN(LEV-1)+THETA_halo_SIZE
      END DO

!  Set pointers for secondary fields in D1
      JEXNER_THETA_LEVELS(1)= SI(406,Sect_No,im_index)
      JP(1)                 = SI(407,Sect_No,im_index)
      JP_THETA_LEVELS(1)    = SI(408,Sect_No,im_index)
      JPSTAR                = SI(409,Sect_No,im_index)
      JSW_INCS(0)           = SI(410,Sect_No,im_index)
      JLW_INCS(0)           = SI(411,Sect_No,im_index)

      DO LEV=1,MODEL_LEVELS+1
        JSW_INCS(LEV)=JSW_INCS(LEV-1)+THETA_FIELD_SIZE
      END DO
      DO LEV=1,MODEL_LEVELS
        JLW_INCS(LEV)=JLW_INCS(LEV-1)+THETA_FIELD_SIZE
      END DO

! Direct PAR flux for STOCHEM
      JDIRPAR               = SI(460,Sect_no,im_index)

      DO LEV=2,MODEL_LEVELS
        JEXNER_THETA_LEVELS(LEV)=JEXNER_THETA_LEVELS(LEV-1)+            &
     &                             theta_off_size
        JP_THETA_LEVELS(LEV)    =JP_THETA_LEVELS(LEV-1)+theta_off_size
      END DO
      DO LEV=2,MODEL_LEVELS+1
        JP(LEV)                 =JP(LEV-1)+theta_off_size
      END DO

!  Set pointers for ancillary fields in D1 from STASH
!     Soil fields
      JSMCL(1)            = SI(  9,Sect_No,im_index)
      J_DEEP_SOIL_TEMP(1) = SI( 20,Sect_No,im_index)
      JVOL_SMC_WILT       = SI( 40,Sect_No,im_index)
      JVOL_SMC_CRIT       = SI( 41,Sect_No,im_index)
      JVOL_SMC_SAT        = SI( 43,Sect_No,im_index)
      JSAT_SOIL_COND      = SI( 44,Sect_No,im_index)
      JTHERM_CAP          = SI( 46,Sect_No,im_index)
      JTHERM_COND         = SI( 47,Sect_No,im_index)
      JSAT_SOILW_SUCTION  = SI( 48,Sect_No,im_index)
      JCLAPP_HORN         = SI(207,Sect_No,im_index)
      JSTHU(1)            = SI(214,Sect_No,im_index)
      JSTHF(1)            = SI(215,Sect_No,im_index)

      DO LEV=2,ST_LEVELS
        J_DEEP_SOIL_TEMP(LEV)=J_DEEP_SOIL_TEMP(LEV-1)+LAND_FIELD
      ENDDO

      DO LEV=2,SM_LEVELS
        JSMCL(LEV)=JSMCL(LEV-1)+LAND_FIELD
        JSTHU(LEV)=JSTHU(LEV-1)+LAND_FIELD
        JSTHF(LEV)=JSTHF(LEV-1)+LAND_FIELD
      END DO

!     Vegetation Fields
      JZ0          = SI( 26,Sect_No,im_index) ! roughness length
      JVEG_FRAC    = SI( 50,Sect_No,im_index) ! vegetation fraction
      JROOT_DEPTH  = SI( 51,Sect_No,im_index) ! root depth
      JSFA         = SI( 52,Sect_No,im_index) ! snow free albedo
      JMDSA        = SI( 53,Sect_No,im_index) ! deep snow albedo
      JSURF_RESIST = SI( 54,Sect_No,im_index) ! surface resistance
      JSURF_CAP    = SI( 55,Sect_No,im_index) ! surface capacity
      JINFILT      = SI( 56,Sect_No,im_index) ! infiltration factor
      JLAI         = SI(208,Sect_No,im_index) ! leaf area index
      JCANHT       = SI(209,Sect_No,im_index) ! canopy height
      JGS          = SI(213,Sect_No,im_index) ! stomatal conductance

! Orography fields
      JOROG_SIL      = SI(17,Sect_No,im_index)   ! Silhouette area
      JOROG_HO2      = SI(18,Sect_No,im_index)   ! Peak to trough ht.
      JOROG_SD       = SI(34,Sect_No,im_index)
      JOROG_GRAD_X   = SI( 5,Sect_No,im_index)
      JOROG_GRAD_Y   = SI( 6,Sect_No,im_index)
      JOROG_GRAD_XX  = SI(35,Sect_No,im_index)
      JOROG_GRAD_XY  = SI(36,Sect_No,im_index)
      JOROG_GRAD_YY  = SI(37,Sect_No,im_index)

! Sea/Sea Ice fields
      JU_SEA         = SI( 28,Sect_No,im_index)
      JV_SEA         = SI( 29,Sect_No,im_index)
      JICE_FRACTION  = SI( 31,Sect_No,im_index)
      JICE_THICKNESS = SI( 32,Sect_No,im_index)
      JTI            = SI( 49,Sect_No,im_index)
      JICE_FRACT_CAT = SI(413,Sect_No,im_index)
      JICE_THICK_CAT = SI(414,Sect_No,im_index)
      JTI_CAT        = SI(415,Sect_No,im_index)
      JU_0_P         = SI(269,Sect_No,im_index)
      JV_0_P         = SI(270,Sect_No,im_index)

! Snow fields
      JSNODEP        = SI( 23,Sect_No,im_index) ! Snow depth over land
      JSNODEP_SEA    = SI( 95,Sect_No,im_index) ! Snow depth on sea ice
      JSNODEP_SEA_CAT= SI(416,Sect_No,im_index) ! Snow depth on ice cats
      JSNSOOT        = SI(221,Sect_No,im_index) ! Snow soot content
      JCATCH_SNOW    = SI(241,Sect_No,im_index)
      JSNOW_GRND     = SI(242,Sect_No,im_index)

! Ozone
      JOZONE(1)     = SI(60,Sect_No,im_index)
! Check for zonal ozone and calculate pointers accordingly
      LEXPAND_OZONE=.FALSE.
      IF (A_LOOKUP(LBNPT,PPINDEX(60,im_index)) == 1) THEN
        LEXPAND_OZONE = .TRUE.
      ENDIF
      DO LEV=2,OZONE_LEVELS
        IF(LEXPAND_OZONE) THEN
          JOZONE(LEV)=JOZONE(LEV-1)+ROWS
        ELSE
          JOZONE(LEV)=JOZONE(LEV-1)+THETA_FIELD_SIZE
        END IF
      END DO

!! Tropopause-based Ozone
      IF (tpps_ozone_levels >  0) THEN
        JTPPSOZONE(1)     = SI(341,Sect_No,im_index)

        !Check for zonal tpps_ozone and calculate pointers accordingly
        LEXPAND_TPPS_OZONE=.FALSE.
        IF (A_LOOKUP(LBNPT,PPINDEX(341,im_index)) == 1) THEN
          LEXPAND_TPPS_OZONE = .TRUE.
        ENDIF
      END IF

      DO LEV=2,TPPS_OZONE_LEVELS
        IF(LEXPAND_TPPS_OZONE) THEN
          JTPPSOZONE(LEV)=JTPPSOZONE(LEV-1)+ROWS
        ELSE
          JTPPSOZONE(LEV)=JTPPSOZONE(LEV-1)+THETA_FIELD_SIZE
        END IF
      END DO

! Add prognostic ozone tracer and cariolle parameters to section 0
! 
        JOZONE_TRACER = SI(480,sect_no,im_index)
        JO3_PROD_LOSS = SI(481,sect_no,im_index)
        JO3_P_L_VMR   = SI(482,sect_no,im_index)
        JO3_VMR       = SI(483,sect_no,im_index)
        JO3_P_L_TEMP  = SI(484,sect_no,im_index)
        JO3_TEMP      = SI(485,sect_no,im_index)
        JO3_P_L_COLO3 = SI(486,sect_no,im_index)
        JO3_COLO3     = SI(487,sect_no,im_index)

        DO  LEV = 2, model_levels
          JOZONE_TRACER(LEV) = JOZONE_TRACER(LEV-1) + THETA_OFF_SIZE
          JO3_PROD_LOSS(LEV) = JO3_PROD_LOSS(LEV-1) + ROWS
          JO3_P_L_VMR(LEV)   = JO3_P_L_VMR(LEV-1) + ROWS
          JO3_VMR(LEV)       = JO3_VMR(LEV-1) + ROWS
          JO3_P_L_TEMP(LEV)  = JO3_P_L_TEMP(LEV-1) + ROWS 
          JO3_TEMP(LEV)      = JO3_TEMP(LEV-1) + ROWS
          JO3_P_L_COLO3(LEV) = JO3_P_L_COLO3(LEV-1) + ROWS
          JO3_COLO3(LEV)     = JO3_COLO3(LEV-1) + ROWS

        END DO


! STOCHEM fields
! Add stochem fields
      JCH4_STOCH(1)=SI(99,Sect_No,im_index)
      JO3_STOCH(1)=SI(100,Sect_No,im_index)

      DO LEV=2,MODEL_LEVELS
        JCH4_STOCH(LEV)=JCH4_STOCH(LEV-1)+THETA_FIELD_SIZE
        JO3_STOCH(LEV)=JO3_STOCH(LEV-1)+THETA_FIELD_SIZE
      ENDDO

! Add sources and aerosol ancillaries

      JMURK_SOURCE(1) = SI(57,Sect_No,im_index) ! Murk source
      JSO2_EM   = SI(58,Sect_No,im_index)   ! Sulphur dioxide emiss.
      JDMS_EM   = SI(59,Sect_No,im_index)   ! Dimethyl sulphide emiss.
      JMURK(1)  = SI(90,Sect_No,im_index)   ! Murk concentration

! Add for Sulphur Cycle
      JSO2(1)       =SI(101,Sect_No,im_index) !Sulphur dioxide gas
      JDMS(1)       =SI(102,Sect_No,im_index) !Dimethyl sulphide gas
      JSO4_AITKEN(1)=SI(103,Sect_No,im_index) !Aitken mode SO4 aerosol
      JSO4_ACCU(1)  =SI(104,Sect_No,im_index) !Accumulation mode SO4 aer
      JSO4_DISS(1)  =SI(105,Sect_No,im_index) !Dissolved SO4 aerosol
      JH2O2(1)      =SI(106,Sect_No,im_index) !Hydrogen peroxide mmr
      JNH3(1)       =SI(107,Sect_No,im_index)  !Ammonia gas
      JSOOT_NEW(1)  =SI(108,Sect_No,im_index)  !Fresh soot
      JSOOT_AGD(1)  =SI(109,Sect_No,im_index)  !Aged soot
      JSOOT_CLD(1)  =SI(110,Sect_No,im_index)  !Soot in cloud
      JBMASS_NEW(1) =SI(111,Sect_No,im_index)  !Fresh biomass smoke
      JBMASS_AGD(1) =SI(112,Sect_No,im_index)  !Aged biomass smoke
      JBMASS_CLD(1) =SI(113,Sect_No,im_index)  !Cloud biomass smoke
      JOCFF_NEW(1)  =SI(114,Sect_No,im_index)  !Fresh ocff
      JOCFF_AGD(1)  =SI(115,Sect_No,im_index)  !Aged socff
      JOCFF_CLD(1)  =SI(116,Sect_No,im_index)  !Ocff in cloud
      JSO2_NATEM(1) =SI(121,Sect_No,im_index)  !Natural SO2 emissions
      JOH(1)        =SI(122,Sect_No,im_index)  !OH 3_D ancillary
      JHO2(1)       =SI(123,Sect_No,im_index)  !HO2 3_D ancillary
      JH2O2_LIMIT(1)=SI(124,Sect_No,im_index)  !H2O2 LIMIT 3_D ancillary
      JO3_CHEM(1)   =SI(125,Sect_No,im_index)  !O3 for chemistry 3_D anc
      JSO2_HILEM    =SI(126,Sect_No,im_index)  !High level SO2 emissions
      JNH3_EM       =SI(127,Sect_No,im_index)  !Ammonia surface emiss
      JSOOT_EM      =SI(128,Sect_No,im_index)  !Fresh soot surf emiss
      JSOOT_HILEM   =SI(129,Sect_No,im_index)  !Fresh soot high emiss
      JBMASS_EM     =SI(130,Sect_No,im_index)  !Fresh bmass surf emiss
      JBMASS_HILEM  =SI(131,Sect_No,im_index)  !Elevated bmass emiss
      JDMS_CONC     =SI(132,Sect_No,im_index)  !DMS conc in seawater
      JDMS_OFLUX    =SI(133,Sect_No,im_index)  !DMS flux from ocean
      JOCFF_EM      =SI(134,Sect_No,im_index)  !Fresh OCFF surf emiss
      JOCFF_HILEM   =SI(135,Sect_No,im_index)  !Fresh OCFF high emiss

! Aerosol climatologies
      JARCLBIOG_BG  =SI(351,Sect_No,im_index)  ! Biogenic aerosol climatology
      JARCLBIOM_FR  =SI(352,Sect_No,im_index)  ! Biomass burning (fresh) aerosol clim
      JARCLBIOM_AG  =SI(353,Sect_No,im_index)  ! Biomass burning (aged) aerosol clim
      JARCLBIOM_IC  =SI(354,Sect_No,im_index)  ! Biomass burning (in-cloud) aerosol clim
      JARCLBLCK_FR  =SI(355,Sect_No,im_index)  ! Black carbon (fresh) aerosol clim
      JARCLBLCK_AG  =SI(356,Sect_No,im_index)  ! Black carbon (aged) aerosol clim
      JARCLSSLT_FI  =SI(357,Sect_No,im_index)  ! Sea salt (film mode) aerosol clim 
      JARCLSSLT_JT  =SI(358,Sect_No,im_index)  ! Sea salt (jet mode) aerosol clim
      JARCLSULP_AC  =SI(359,Sect_No,im_index)  ! Sulphate (accumulation mode) aero clim
      JARCLSULP_AK  =SI(360,Sect_No,im_index)  ! Sulphate (Aitken mode) aerosol clim 
      JARCLSULP_DI  =SI(361,Sect_No,im_index)  ! Sulphate (dissolved) aerosol clim
      JARCLDUST_B1  =SI(362,Sect_No,im_index)  ! Dust (bin 1) aerosol climatology 
      JARCLDUST_B2  =SI(363,Sect_No,im_index)  ! Dust (bin 2) aerosol climatology 
      JARCLDUST_B3  =SI(364,Sect_No,im_index)  ! Dust (bin 3) aerosol climatology 
      JARCLDUST_B4  =SI(365,Sect_No,im_index)  ! Dust (bin 4) aerosol climatology 
      JARCLDUST_B5  =SI(366,Sect_No,im_index)  ! Dust (bin 5) aerosol climatology 
      JARCLDUST_B6  =SI(367,Sect_No,im_index)  ! Dust (bin 6) aerosol climatology 
      JARCLOCFF_FR  =SI(368,Sect_No,im_index)  ! Org carbon fossil fuel (fresh) aero clim
      JARCLOCFF_AG  =SI(369,Sect_No,im_index)  ! Org carbon fossil fuel (aged) aero clim
      JARCLOCFF_IC  =SI(370,Sect_No,im_index)  ! Org carbon fossil fuel (in-cloud) aero clim
      JARCLDLTA_DL  =SI(371,Sect_No,im_index)  ! Delta aerosol climatology

      DO LEV=2,MODEL_LEVELS
        JARCLBIOG_BG(LEV) = JARCLBIOG_BG(LEV-1)+THETA_FIELD_SIZE
        JARCLBIOM_FR(LEV) = JARCLBIOM_FR(LEV-1)+THETA_FIELD_SIZE
        JARCLBIOM_AG(LEV) = JARCLBIOM_AG(LEV-1)+THETA_FIELD_SIZE
        JARCLBIOM_IC(LEV) = JARCLBIOM_IC(LEV-1)+THETA_FIELD_SIZE
        JARCLBLCK_FR(LEV) = JARCLBLCK_FR(LEV-1)+THETA_FIELD_SIZE
        JARCLBLCK_AG(LEV) = JARCLBLCK_AG(LEV-1)+THETA_FIELD_SIZE
        JARCLSSLT_FI(LEV) = JARCLSSLT_FI(LEV-1)+THETA_FIELD_SIZE
        JARCLSSLT_JT(LEV) = JARCLSSLT_JT(LEV-1)+THETA_FIELD_SIZE
        JARCLSULP_AC(LEV) = JARCLSULP_AC(LEV-1)+THETA_FIELD_SIZE
        JARCLSULP_AK(LEV) = JARCLSULP_AK(LEV-1)+THETA_FIELD_SIZE
        JARCLSULP_DI(LEV) = JARCLSULP_DI(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B1(LEV) = JARCLDUST_B1(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B2(LEV) = JARCLDUST_B2(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B3(LEV) = JARCLDUST_B3(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B4(LEV) = JARCLDUST_B4(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B5(LEV) = JARCLDUST_B5(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B6(LEV) = JARCLDUST_B6(LEV-1)+THETA_FIELD_SIZE
        JARCLOCFF_FR(LEV) = JARCLOCFF_FR(LEV-1)+THETA_FIELD_SIZE
        JARCLOCFF_AG(LEV) = JARCLOCFF_AG(LEV-1)+THETA_FIELD_SIZE
        JARCLOCFF_IC(LEV) = JARCLOCFF_IC(LEV-1)+THETA_FIELD_SIZE
        JARCLDLTA_DL(LEV) = JARCLDLTA_DL(LEV-1)+THETA_FIELD_SIZE
      END DO

! Mineral dust scheme

      JSOIL_CLAY   =SI(418,Sect_No,im_index)  ! soil clay fraction
      JSOIL_SILT   =SI(419,Sect_No,im_index)  ! soil silt fraction
      JSOIL_SAND   =SI(420,Sect_No,im_index)  ! soil sand fraction

      JDUST_MREL1=SI(421,Sect_No,im_index) !relative soil mass in div1
      JDUST_MREL2=SI(422,Sect_No,im_index) !relative soil mass in div2
      JDUST_MREL3=SI(423,Sect_No,im_index) !relative soil mass in div3
      JDUST_MREL4=SI(424,Sect_No,im_index) !relative soil mass in div4
      JDUST_MREL5=SI(425,Sect_No,im_index) !relative soil mass in div5
      JDUST_MREL6=SI(426,Sect_No,im_index) !relative soil mass in div6


      JDUST_DIV1(1)=SI(431,Sect_No,im_index)  ! dust mmr, division 1
      JDUST_DIV2(1)=SI(432,Sect_No,im_index)  ! dust mmr, division 2
      JDUST_DIV3(1)=SI(433,Sect_No,im_index)  ! dust mmr, division 3
      JDUST_DIV4(1)=SI(434,Sect_No,im_index)  ! dust mmr, division 4
      JDUST_DIV5(1)=SI(435,Sect_No,im_index)  ! dust mmr, division 5
      JDUST_DIV6(1)=SI(436,Sect_No,im_index)  ! dust mmr, division 6

      DO LEV = 2,MODEL_LEVELS
       JDUST_DIV1(LEV)=JDUST_DIV1(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV2(LEV)=JDUST_DIV2(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV3(LEV)=JDUST_DIV3(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV4(LEV)=JDUST_DIV4(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV5(LEV)=JDUST_DIV5(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV6(LEV)=JDUST_DIV6(LEV-1)+THETA_OFF_SIZE
      ENDDO

! HadCM2 sulphate loading patterns
      JHadCM2_SO4(1)=SI(160,Sect_No,im_index)
      DO LEV=2, NSULPAT
        JHadCM2_SO4(LEV)=JHadCM2_SO4(LEV-1)+THETA_FIELD_SIZE
      ENDDO

! Add for Carbon cycle
      J_CO2FLUX = SI(250,Sect_No,im_index)
      J_CO2_EMITS  = SI(251,Sect_No,im_index)
      JCO2(1)      = SI(252,Sect_No,im_index)
      DO LEV=2,MODEL_LEVELS
        JMURK_SOURCE(LEV) = JMURK_SOURCE(LEV-1)+THETA_FIELD_SIZE
        JMURK(LEV) = JMURK(LEV-1)+THETA_OFF_SIZE

! For Sulphur Cycle variables
        JSO2(LEV)=JSO2(LEV-1)+THETA_OFF_SIZE
        JDMS(LEV)=JDMS(LEV-1)+THETA_OFF_SIZE
        JSO4_AITKEN(LEV)=JSO4_AITKEN(LEV-1)+THETA_OFF_SIZE
        JSO4_ACCU(LEV)=JSO4_ACCU(LEV-1)+THETA_OFF_SIZE
        JSO4_DISS(LEV)=JSO4_DISS(LEV-1)+THETA_OFF_SIZE
        JH2O2(LEV)=JH2O2(LEV-1)+THETA_OFF_SIZE
        JSO2_NATEM(LEV)=JSO2_NATEM(LEV-1)+THETA_FIELD_SIZE
        JOH(LEV) = JOH(LEV-1)+THETA_FIELD_SIZE
        JHO2(LEV) = JHO2(LEV-1)+THETA_FIELD_SIZE
        JNH3(LEV)      = JNH3(LEV-1)+THETA_OFF_SIZE
        JSOOT_NEW(LEV) = JSOOT_NEW(LEV-1)+THETA_OFF_SIZE
        JSOOT_AGD(LEV) = JSOOT_AGD(LEV-1)+THETA_OFF_SIZE
        JSOOT_CLD(LEV) = JSOOT_CLD(LEV-1)+THETA_OFF_SIZE
        JBMASS_NEW(LEV) = JBMASS_NEW(LEV-1)+THETA_OFF_SIZE
        JBMASS_AGD(LEV) = JBMASS_AGD(LEV-1)+THETA_OFF_SIZE
        JBMASS_CLD(LEV) = JBMASS_CLD(LEV-1)+THETA_OFF_SIZE
        JOCFF_NEW(LEV) = JOCFF_NEW(LEV-1)+THETA_OFF_SIZE
        JOCFF_AGD(LEV) = JOCFF_AGD(LEV-1)+THETA_OFF_SIZE
        JOCFF_CLD(LEV) = JOCFF_CLD(LEV-1)+THETA_OFF_SIZE
        JH2O2_LIMIT(LEV)=JH2O2_LIMIT(LEV-1)+THETA_FIELD_SIZE
        JO3_CHEM(LEV)=JO3_CHEM(LEV-1)+THETA_FIELD_SIZE

! For Carbon Cycle variables
        JCO2(LEV)=JCO2(LEV-1)+THETA_OFF_SIZE
      END DO

! Tracer fluxes - kdcorbin, 06/10
      tflux_sect = 3
      JTRACER_FLUX1 = SI(100,tflux_sect,im_index)
      JTRACER_FLUX2 = SI(101,tflux_sect,im_index)
      JTRACER_FLUX3 = SI(102,tflux_sect,im_index)
      JTRACER_FLUX4 = SI(103,tflux_sect,im_index)
      JTRACER_FLUX5 = SI(104,tflux_sect,im_index)
      JTRACER_FLUX6 = SI(105,tflux_sect,im_index)
      JTRACER_FLUX7 = SI(106,tflux_sect,im_index)
      JTRACER_FLUX8 = SI(107,tflux_sect,im_index)
      JTRACER_FLUX9 = SI(108,tflux_sect,im_index)
      JTRACER_FLUX10= SI(109,tflux_sect,im_index)
      JTRACER_FLUX11= SI(110,tflux_sect,im_index)
      JTRACER_FLUX12= SI(111,tflux_sect,im_index)
      JTRACER_FLUX13= SI(112,tflux_sect,im_index)
      JTRACER_FLUX14= SI(113,tflux_sect,im_index)
      JTRACER_FLUX15= SI(114,tflux_sect,im_index)
      JTRACER_FLUX16= SI(115,tflux_sect,im_index)
      JTRACER_FLUX17= SI(116,tflux_sect,im_index)
      JTRACER_FLUX18= SI(117,tflux_sect,im_index)
      JTRACER_FLUX19= SI(118,tflux_sect,im_index)
      JTRACER_FLUX20= SI(119,tflux_sect,im_index)

!jhan(re-assign atm_pointers
! Lest Sept 2012 - changed items from 300s (user ancils) to 800s
      !--- allow for 6 layers
      JTSOIL_TILE(1)   = SI(801,Sect_No,im_index)
      JTSOIL_TILE(2)   = SI(802,Sect_No,im_index)
      JTSOIL_TILE(3)   = SI(803,Sect_No,im_index)
      JTSOIL_TILE(4)   = SI(804,Sect_No,im_index)
      JTSOIL_TILE(5)   = SI(805,Sect_No,im_index)
      JTSOIL_TILE(6)   = SI(806,Sect_No,im_index)
      JSMCL_TILE(1)    = SI(807,Sect_No,im_index)
      JSMCL_TILE(2)    = SI(808,Sect_No,im_index)
      JSMCL_TILE(3)    = SI(809,Sect_No,im_index)
      JSMCL_TILE(4)    = SI(810,Sect_No,im_index)
      JSMCL_TILE(5)    = SI(811,Sect_No,im_index)
      JSMCL_TILE(6)    = SI(812,Sect_No,im_index)
      JSTHF_TILE(1)    = SI(813,Sect_No,im_index)
      JSTHF_TILE(2)    = SI(814,Sect_No,im_index)
      JSTHF_TILE(3)    = SI(815,Sect_No,im_index)
      JSTHF_TILE(4)    = SI(816,Sect_No,im_index)
      JSTHF_TILE(5)    = SI(817,Sect_No,im_index)
      JSTHF_TILE(6)    = SI(818,Sect_No,im_index)
      JSNOW_DEPTH3L(1) = SI(819,Sect_No,im_index)
      JSNOW_DEPTH3L(2) = SI(820,Sect_No,im_index)
      JSNOW_DEPTH3L(3) = SI(821,Sect_No,im_index)
      JSNOW_MASS3L(1)  = SI(822,Sect_No,im_index)
      JSNOW_MASS3L(2)  = SI(823,Sect_No,im_index)
      JSNOW_MASS3L(3)  = SI(824,Sect_No,im_index)
      JSNOW_TMP3L(1)   = SI(825,Sect_No,im_index)
      JSNOW_TMP3L(2)   = SI(826,Sect_No,im_index)
      JSNOW_TMP3L(3)   = SI(827,Sect_No,im_index)
      JSNOW_RHO3L(1)   = SI(828,Sect_No,im_index)
      JSNOW_RHO3L(2)   = SI(829,Sect_No,im_index)
      JSNOW_RHO3L(3)   = SI(830,Sect_No,im_index)
      JSNOW_RHO1L      = SI(831,Sect_No,im_index)
      JSNOW_AGE        = SI(832,Sect_No,im_index)
      JSNOW_FLG3l      = SI(833,Sect_No,im_index)

! Lestevens Sept 2012 - START CASA
! sect_no equals 0
      JCPOOL_TILE(1)  =SI(851,Sect_No,im_index)  !
      JCPOOL_TILE(2)  =SI(852,Sect_No,im_index)  !
      JCPOOL_TILE(3)  =SI(853,Sect_No,im_index)  !
      JCPOOL_TILE(4)  =SI(854,Sect_No,im_index)  !
      JCPOOL_TILE(5)  =SI(855,Sect_No,im_index)  !
      JCPOOL_TILE(6)  =SI(856,Sect_No,im_index)  !
      JCPOOL_TILE(7)  =SI(857,Sect_No,im_index)  !
      JCPOOL_TILE(8)  =SI(858,Sect_No,im_index)  !
      JCPOOL_TILE(9)  =SI(859,Sect_No,im_index)  !
      JCPOOL_TILE(10) =SI(860,Sect_No,im_index)  !
      JNPOOL_TILE(1)  =SI(861,Sect_No,im_index)  !
      JNPOOL_TILE(2)  =SI(862,Sect_No,im_index)  !
      JNPOOL_TILE(3)  =SI(863,Sect_No,im_index)  !
      JNPOOL_TILE(4)  =SI(864,Sect_No,im_index)  !
      JNPOOL_TILE(5)  =SI(865,Sect_No,im_index)  !
      JNPOOL_TILE(6)  =SI(866,Sect_No,im_index)  !
      JNPOOL_TILE(7)  =SI(867,Sect_No,im_index)  !
      JNPOOL_TILE(8)  =SI(868,Sect_No,im_index)  !
      JNPOOL_TILE(9)  =SI(869,Sect_No,im_index)  !
      JNPOOL_TILE(10) =SI(870,Sect_No,im_index)  !
      JPPOOL_TILE(1)  =SI(871,Sect_No,im_index)  !
      JPPOOL_TILE(2)  =SI(872,Sect_No,im_index)  !
      JPPOOL_TILE(3)  =SI(873,Sect_No,im_index)  !
      JPPOOL_TILE(4)  =SI(874,Sect_No,im_index)  !
      JPPOOL_TILE(5)  =SI(875,Sect_No,im_index)  !
      JPPOOL_TILE(6)  =SI(876,Sect_No,im_index)  !
      JPPOOL_TILE(7)  =SI(877,Sect_No,im_index)  !
      JPPOOL_TILE(8)  =SI(878,Sect_No,im_index)  !
      JPPOOL_TILE(9)  =SI(879,Sect_No,im_index)  !
      JPPOOL_TILE(10) =SI(880,Sect_No,im_index)  !
      JPPOOL_TILE(11) =SI(881,Sect_No,im_index)  !
      JPPOOL_TILE(12) =SI(882,Sect_No,im_index)  !
      JSOIL_ORDER     =SI(883,Sect_No,im_index)  !
      JNIDEP          =SI(884,Sect_No,im_index)  !
      JNIFIX          =SI(885,Sect_No,im_index)  !
      JPDUST          =SI(887,Sect_No,im_index)  !
      JPWEA           =SI(888,Sect_No,im_index)  !
      JGLAI           =SI(893,Sect_No,im_index)  !
      JPHENPHASE      =SI(894,Sect_No,im_index)  !

!     JNFERT          =SI(886,Sect_No,im_index)  !
!     JPFERT          =SI(889,Sect_No,im_index)  !
!     JCBAL_SUM       =SI(890,Sect_No,im_index)  !
!     JNBAL_SUM       =SI(891,Sect_No,im_index)  !
!     JPBAL_SUM       =SI(892,Sect_No,im_index)  !
! Lestevens Sept 2012 - END CASA

      JCABLE_LAI(1)   = SI(895,Sect_No,im_index)   

! Add user ancillaries

      JUSER_ANC1  = SI(301,Sect_No,im_index)
      JUSER_ANC2  = SI(302,Sect_No,im_index)
      JUSER_ANC3  = SI(303,Sect_No,im_index)
      JUSER_ANC4  = SI(304,Sect_No,im_index)
      JUSER_ANC5  = SI(305,Sect_No,im_index)
      JUSER_ANC6  = SI(306,Sect_No,im_index)
      JUSER_ANC7  = SI(307,Sect_No,im_index)
      JUSER_ANC8  = SI(308,Sect_No,im_index)
      JUSER_ANC9  = SI(309,Sect_No,im_index)
      JUSER_ANC10 = SI(310,Sect_No,im_index)
      JUSER_ANC11 = SI(311,Sect_No,im_index)
      JUSER_ANC12 = SI(312,Sect_No,im_index)
      JUSER_ANC13 = SI(313,Sect_No,im_index)
      JUSER_ANC14 = SI(314,Sect_No,im_index)
      JUSER_ANC15 = SI(315,Sect_No,im_index)
      JUSER_ANC16 = SI(316,Sect_No,im_index)
      JUSER_ANC17 = SI(317,Sect_No,im_index)
      JUSER_ANC18 = SI(318,Sect_No,im_index)
      JUSER_ANC19 = SI(319,Sect_No,im_index)
      JUSER_ANC20 = SI(320,Sect_No,im_index)
      JUSER_MULT1(1)  = SI(321,Sect_No,im_index)
      JUSER_MULT2(1)  = SI(322,Sect_No,im_index)
      JUSER_MULT3(1)  = SI(323,Sect_No,im_index)
      JUSER_MULT4(1)  = SI(324,Sect_No,im_index)
      JUSER_MULT5(1)  = SI(325,Sect_No,im_index)
      JUSER_MULT6(1)  = SI(326,Sect_No,im_index)
      JUSER_MULT7(1)  = SI(327,Sect_No,im_index)
      JUSER_MULT8(1)  = SI(328,Sect_No,im_index)
      JUSER_MULT9(1)  = SI(329,Sect_No,im_index)
      JUSER_MULT10(1) = SI(330,Sect_No,im_index)
      JUSER_MULT11(1) = SI(331,Sect_No,im_index)
      JUSER_MULT12(1) = SI(332,Sect_No,im_index)
      JUSER_MULT13(1) = SI(333,Sect_No,im_index)
      JUSER_MULT14(1) = SI(334,Sect_No,im_index)
      JUSER_MULT15(1) = SI(335,Sect_No,im_index)
      JUSER_MULT16(1) = SI(336,Sect_No,im_index)
      JUSER_MULT17(1) = SI(337,Sect_No,im_index)
      JUSER_MULT18(1) = SI(338,Sect_No,im_index)
      JUSER_MULT19(1) = SI(339,Sect_No,im_index)
      JUSER_MULT20(1) = SI(340,Sect_No,im_index)

! Set for multi-level user ancillaries
      DO LEV=2,MODEL_LEVELS
        JUSER_MULT1(LEV)  = JUSER_MULT1(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT2(LEV)  = JUSER_MULT2(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT3(LEV)  = JUSER_MULT3(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT4(LEV)  = JUSER_MULT4(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT5(LEV)  = JUSER_MULT5(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT6(LEV)  = JUSER_MULT6(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT7(LEV)  = JUSER_MULT7(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT8(LEV)  = JUSER_MULT8(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT9(LEV)  = JUSER_MULT9(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT10(LEV) = JUSER_MULT10(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT11(LEV) = JUSER_MULT11(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT12(LEV) = JUSER_MULT12(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT13(LEV) = JUSER_MULT13(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT14(LEV) = JUSER_MULT14(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT15(LEV) = JUSER_MULT15(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT16(LEV) = JUSER_MULT16(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT17(LEV) = JUSER_MULT17(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT18(LEV) = JUSER_MULT18(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT19(LEV) = JUSER_MULT19(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT20(LEV) = JUSER_MULT20(LEV-1)+THETA_FIELD_SIZE
      END DO

! Tiled vegetation and triffid
      JFRAC_TYP     = SI(216,Sect_No,im_index) ! surface type fractions
      JFRAC_CON1    = SI(442,Sect_No,im_index) ! surface type fractions
      JFRAC_CON2    = SI(443,Sect_No,im_index) ! surface type fractions
      JFRAC_CON3    = SI(444,Sect_No,im_index) ! surface type fractions
      JFRAC_CON4    = SI(445,Sect_No,im_index) ! surface type fractions
      JFRAC_CON5    = SI(446,Sect_No,im_index) ! surface type fractions
      JFRAC_CON6    = SI(447,Sect_No,im_index) ! surface type fractions
      JFRAC_CON7    = SI(448,Sect_No,im_index) ! surface type fractions
      JFRAC_CON8    = SI(449,Sect_No,im_index) ! surface type fractions
      JFRAC_CON9    = SI(450,Sect_No,im_index) ! surface type fractions
      JLAI_PFT      = SI(217,Sect_No,im_index) ! leaf area index of PFTs
      JCANHT_PFT    = SI(218,Sect_No,im_index) ! canopy height of PFTs
      JDISTURB      = SI(219,Sect_No,im_index) ! Veg disturbed fraction
      JSOIL_ALB     = SI(220,Sect_No,im_index) ! Snow-free soil albedo
      JSOIL_CARB    = SI(223,Sect_No,im_index) ! Soil carbon content
      JSOIL_CARB1   = SI(466,Sect_No,im_index) ! Soil carbon content DPM
      JSOIL_CARB2   = SI(467,Sect_No,im_index) ! Soil carbon content RPM
      JSOIL_CARB3   = SI(468,Sect_No,im_index) ! Soil carbon content BIO
      JSOIL_CARB4   = SI(469,Sect_No,im_index) ! Soil carbon content HUM
      JNPP_PFT_ACC  = SI(224,Sect_No,im_index) ! Accumulated NPP on PFTs
      JG_LF_PFT_ACC = SI(225,Sect_No,im_index) ! Accumulated leaf
!                                              ! turnover rate on PFTs
      JG_PHLF_PFT_ACC=SI(226,Sect_No,im_index) ! Accumulat. phenological
!                                              ! leaf turnover rate PFTs
      JRSP_W_PFT_ACC= SI(227,Sect_No,im_index) ! Accum. wood resp PFTs
      JRSP_S_ACC    = SI(228,Sect_No,im_index) ! Accumulated soil resp
      JRSP_S_ACC1 = SI(470,Sect_No,im_index)  ! Soil respiration DPM
      JRSP_S_ACC2 = SI(471,Sect_No,im_index)  ! Soil respiration RPM
      JRSP_S_ACC3 = SI(472,Sect_No,im_index)  ! Soil respiration BIO
      JRSP_S_ACC4 = SI(473,Sect_No,im_index)  ! Soil respiration HUM
      JCAN_WATER_TILE=SI(229,Sect_No,im_index) ! Canopy water content
!                                              ! on tiles
      JCATCH_TILE   = SI(230,Sect_No,im_index) ! Canopy capacity on
!                                              ! tiles
      JRGRAIN_TILE  = SI(231,Sect_No,im_index) ! Snow grain size on
!                                              ! tiles
      JTSTAR_TILE   = SI(233,Sect_No,im_index) ! Tiled surface temp
      JZ0_TILE      = SI(234,Sect_No,im_index) ! Tiled surface roughness
! Stash number for snow tile not finalised yet
      JSNODEP_TILE  = SI(240,Sect_No,im_index) ! Tiled snow depth
      JINFIL_TILE   = SI(236,Sect_No,im_index) ! Max tile infilt rate
! Stash codes for DORL, LW_DOWN, SW_TILE not finalised yet
      JDOLR         = SI(239,Sect_No,im_index) ! TOA surface up LW
      JLW_DOWN      = SI(238,Sect_No,im_index) ! Surface down LW
      JSW_TILE      = SI(237,Sect_No,im_index) ! Surface net SW on tiles


! River routing fields
      JRIV_SEQUENCE  = SI(151,Sect_No,im_index) ! River sequence
      JRIV_DIRECTION = SI(152,Sect_No,im_index) ! River Direction
      JRIV_STORAGE   = SI(153,Sect_No,im_index) ! River Water Storage
      JTOT_SURFROFF  = SI(155,Sect_No,im_index) ! Acc. surface runoff
      JTOT_SUBROFF   = SI(156,Sect_No,im_index) ! Acc. sub-surf runoff
! Set pointer for inland basin outflow
      JRIV_INLANDATM    = SI(511,Sect_No,im_index)
       !Inland basin outflow



! required for energy correction
      JNET_FLUX=SI(222,Sect_No,im_index)    ! store for energy flux
      JNET_MFLUX=SI(235,Sect_No,im_index)   ! store for moisture flux

!Fields carried forward from previous version.
      JTSTAR_ANOM    = SI(39,Sect_No,im_index)




! Set all pointers referencing Lateral Boundary Conditions

      Sect_No=31  ! LBC section

      JOROG_LBC     = SI(1,Sect_No,im_index)
      JU_LBC        = SI(2,Sect_No,im_index)
      JV_LBC        = SI(3,Sect_No,im_index)
      JW_LBC        = SI(4,Sect_No,im_index)
      JRHO_LBC      = SI(5,Sect_No,im_index)
      JTHETA_LBC    = SI(6,Sect_No,im_index)
      JQ_LBC        = SI(7,Sect_No,im_index)
      JQCL_LBC      = SI(8,Sect_No,im_index)
      JQCF_LBC      = SI(9,Sect_No,im_index)
      JEXNER_LBC    = SI(10,Sect_No,im_index)
      JU_ADV_LBC    = SI(11,Sect_No,im_index)
      JV_ADV_LBC    = SI(12,Sect_No,im_index)
      JW_ADV_LBC    = SI(13,Sect_No,im_index)
      JQCF2_LBC     = SI(14,Sect_No,im_index)
      JQRAIN_LBC    = SI(15,Sect_No,im_index)
      JQGRAUP_LBC   = SI(16,Sect_No,im_index)
      JCF_BULK_LBC  = SI(17,Sect_No,im_index)
      JCF_LIQUID_LBC= SI(18,Sect_No,im_index)
      JCF_FROZEN_LBC= SI(19,Sect_No,im_index)
      JMURK_LBC     = SI(20,Sect_No,im_index)

      JU_LBC_TEND     = SI(257,Sect_No,im_index)
      JV_LBC_TEND     = SI(258,Sect_No,im_index)
      JW_LBC_TEND     = SI(259,Sect_No,im_index)
      JRHO_LBC_TEND   = SI(260,Sect_No,im_index)
      JTHETA_LBC_TEND = SI(261,Sect_No,im_index)
      JQ_LBC_TEND     = SI(262,Sect_No,im_index)
      JQCL_LBC_TEND   = SI(263,Sect_No,im_index)
      JQCF_LBC_TEND   = SI(264,Sect_No,im_index)
      JEXNER_LBC_TEND = SI(265,Sect_No,im_index)
      JU_ADV_LBC_TEND = SI(266,Sect_No,im_index)
      JV_ADV_LBC_TEND = SI(267,Sect_No,im_index)
      JW_ADV_LBC_TEND = SI(268,Sect_No,im_index)
      JQCF2_LBC_TEND   = SI(269,Sect_No,im_index)
      JQRAIN_LBC_TEND  = SI(270,Sect_No,im_index)
      JQGRAUP_LBC_TEND = SI(271,Sect_No,im_index)
      JCF_BULK_LBC_TEND  = SI(272,Sect_No,im_index)
      JCF_LIQUID_LBC_TEND= SI(273,Sect_No,im_index)
      JCF_FROZEN_LBC_TEND= SI(274,Sect_No,im_index)
      JMURK_LBC_TEND   = SI(275,Sect_No,im_index)

      IF (TR_VARS  >   0) THEN
        ivar=1
        DO jvar= A_TRACER_FIRST,A_TRACER_LAST
          IF (SI(jvar,33,im_index)  /=  1) THEN
            JTRACER_LBC(ivar)=                                          &
     &        SI(100+jvar-A_TRACER_FIRST,Sect_no,im_index)
            JTRACER_LBC_TEND(ivar)=                                     &
     &        SI(356+jvar-A_TRACER_FIRST,Sect_no,im_index)
            ivar=ivar+1
          ENDIF
        ENDDO ! j
      ELSE
        JTRACER_LBC(1) = 1  ! ensure sensible address even if no tracers
        JTRACER_LBC_TEND(1) = 1
      ENDIF ! IF (TR_VARS  >   0)

! Tracer prognostics are now in section 33, not section 0
      sect_no = 33   ! tracers section
      JVAR=0         ! JVAR+1 is the current tracer to be found
      IF (TR_VARS >  0) THEN
        DO IVAR=A_TRACER_FIRST,A_TRACER_LAST
          IF(SI(IVAR,sect_no,im_index) /= 1) THEN ! tracer in use
            JVAR=JVAR+1
            JTRACER(1,JVAR) = SI(IVAR,sect_no,im_index)
            DO LEV=2,TR_LEVELS
              JTRACER(LEV,JVAR)=JTRACER(LEV-1,JVAR)+THETA_OFF_SIZE
            END DO
            A_TR_INDEX(IVAR-A_TRACER_FIRST+1)=JVAR
          END IF
        END DO
      ELSE
        JTRACER(1,1)=1   ! Ensure a sensible address even if no tracers
      ENDIF
      IF(JVAR /= TR_VARS) THEN
        WRITE(6,*) 'STATMPT: TR_VARS and SI are inconsistent'
        WRITE(6,*) 'TR_VARS=',TR_VARS,' .     But, SI implies :',JVAR
        CMESSAGE=  'STATMPT: TR_VARS and SI  inconsistent, see output'
        ICODE=100
        GOTO 9999 ! error return
      END IF
! Oxidant concentrations from UKCA for use in HadGEM sulphur
! cycle (these are in Section 33):
      IF(L_SULPC_ONLINE_OXIDANTS .AND. L_UKCA) THEN 
         JOH_UKCA(1)   =SI(251,sect_no,im_index) 
         JH2O2_UKCA(1) =SI(8,sect_no,im_index)   
         JHO2_UKCA(1)  =SI(252,sect_no,im_index) 
         JO3_UKCA(1)   =SI(1,sect_no,im_index)   
      END IF 

! UKCA tracer prognostics are in section 34.
      sect_no = 34   ! UKCA tracers section
      JVAR=0         ! JVAR+1 is the current tracer to be found
        DO IVAR=A_UKCA_FIRST,A_UKCA_LAST
          IF(SI(IVAR,sect_no,im_index) /= 1) THEN ! tracer in use
            JVAR=JVAR+1
            JTR_UKCA(1,JVAR) = SI(IVAR,sect_no,im_index)
            DO LEV=2,TR_LEVELS
              JTR_UKCA(LEV,JVAR)=JTR_UKCA(LEV-1,JVAR)+THETA_OFF_SIZE
            END DO
            A_UKCA_INDEX(IVAR-A_UKCA_FIRST+1)=JVAR
          END IF
        END DO
      IF (JVAR  ==  0) THEN
        JTR_UKCA(1,1)=1   ! Ensure a sensible address when no tracers
      ENDIF

! Set pointers to level dependent constants for atmosphere.
      JETATHETA      =1              ! eta_theta_levels(0:model_levels)
      JETARHO        =JETATHETA+model_levels+1 ! eta_rho_levels
      JRHCRIT        =JETARHO+model_levels+1   ! rhcrit
      JSOIL_THICKNESS=JRHCRIT+model_levels+1   ! soil level depths
! For height definition z(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      Jzseak_theta   =JSOIL_THICKNESS+model_levels+1 ! zsea (theta levs)
      JCk_theta      =Jzseak_theta   +model_levels+1 ! C    (theta levs)
      Jzseak_rho     =JCk_theta      +model_levels+1 ! zsea (rho levs)
      JCk_rho        =Jzseak_rho     +model_levels+1 ! C    (rho levs)
! Set pointers to Row dependent constants for atmosphere.
      JPHI_INPUT_P       =1
      JPHI_INPUT_V       =JPHI_INPUT_P + global_rows
! Set pointers to Col dependent constants for atmosphere.      
      JLAMBDA_INPUT_P    =1
      JLAMBDA_INPUT_U    =JLAMBDA_INPUT_P + global_row_length

! Set pointers to Row and Col dependent constants for atmosphere.
      JPHI_INPUT_P       =1
      JPHI_INPUT_V       =JPHI_INPUT_P + global_rows
      JLAMBDA_INPUT_P    =1
      JLAMBDA_INPUT_U    =JLAMBDA_INPUT_P + global_row_length

! The following if block should only be required for test purposes
! during the transition of vn5.2. This ensures old values for
! setting blev/bulev/brlev on old style dumps (before 'smooth'
! algorithm introduced and new definitions introduced).
      if(A_LEN2_LEVDEPC  <=  4) then ! ie before 5.2 change
         JCk_rho=jetarho
         Jck_theta=jetatheta
         Jzseak_rho=jetarho
         Jzseak_theta=jetatheta
         write(6,*) 'SETCONA_CTL: WARNING. This dump has not been '//   &
     &   'reconfigured with new 5.2 set of level dependent constants:'
         write(6,*) 'PP headers will revert to 5.0/5.1 style '//        &
     &   'blev/bulev/brlev definitions'
      endif

 9999 CONTINUE ! ERROR GOTO point.
      RETURN
      END SUBROUTINE SET_ATM_POINTERS

