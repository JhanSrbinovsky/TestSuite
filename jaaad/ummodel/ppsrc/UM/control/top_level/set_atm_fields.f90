
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Points atmosphere fields to the appropriate sections of D1
!
! Subroutine Interface: 
SUBROUTINE Set_Atm_Fields ( &
! "argptra.h" contains jpointers
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
! "argsts.h" contains SI (STASH index array); used to check tracers
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
      D1, LD1, ID1 )

       USE atm_fields_mod   ! atmosphere fields
       USE field_length_mod ! field_length function 

IMPLICIT NONE
!
! Description: 
!   Routine to point atmosphere fields to the appropriate sections of D1.
!   After calling this subroutine, the fields can be used directly without
!   referring to D1
!
! Method: 
!   Assumming SET_ATM_POINTERS has been called beforehand, this subroutine
!   points each field to an area of D1 starting at the corresponding
!   "jpointer" (at the first level) and ending at the "jpointer" (at the 
!   last level) plus the size of a single level of that field. 
!
!   Tracers are dealt with differently:   First the number of active tracers
!   is computed so that the correct sections of the corresponding tracer
!   jpointers can be used in pointing the tracer fields to D1.  If no tracers
!   are active then the fields are pointed to a dummy array
!
! Owner:  A. Treshansky
!
! Code description: 
!   Language:  Fortran 90.
!   This code is written to UM programming standards version 8.0.
!

! Subroutine arguments

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
! constants (N_INTERNAL_MODELS_MAX) in "csubmodl.h" are needed by "ccontrol.h"
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
! constants (NUNITS) in "chsunits.h" are needed by "ccontrol.h"
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
! constants (L_3D_CCA) in "ccontrol.h" are needed to determine field sizes below
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

! constants (A_TRACER_FIRST, etc.) in "ctracera.h" are needed for tracers
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

      REAL,    TARGET, INTENT(IN) :: D1(LEN_TOT)
      LOGICAL, TARGET, INTENT(IN) :: LD1(LEN_TOT)
      INTEGER, TARGET, INTENT(IN) :: ID1(LEN_TOT)

! Local variables

      INTEGER :: nTracer ! loop counter over available tracers
      INTEGER :: nActiveTracers ! number of tracers actually being used

! End of header

! 1.0 Start of subroutine code; point fields to D1

!    atmospheric primary variables
     U         => D1(JU(1) : JU(model_levels)+u_off_size )
     V         => D1(JV(1) : JV(model_levels)+v_off_size )
     THETA     => D1(JTHETA(1) : JTHETA(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     Q         => D1(JQ(1) : JQ(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     QCF       => D1(JQCF(1) : JQCF(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     TSTAR     => D1(JTSTAR : JTSTAR+field_length(theta_points,no_halo,1) )
     LAND      => D1(JLAND : JLAND+field_length(theta_points,no_halo,1) )
     OROGRAPHY => D1(JOROG : JOROG+field_length(theta_points,no_halo,1) )
     W         => D1(JW(0) : JW(model_levels)+theta_off_size)
     RHO       => D1(JRHO(1) : JRHO(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     QCL       => D1(JQCL(1) : JQCL(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     QCF2      => D1(JQCF2(1) : JQCF2(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     QRAIN     => D1(JQRAIN(1) : JQRAIN(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     QGRAUP    => D1(JQGRAUP(1) : JQGRAUP(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     EXNER_RHO_LEVELS => D1(JEXNER_RHO_LEVELS(1) : JEXNER_RHO_LEVELS(1)+ &
      field_length(theta_points,single_halo,model_levels+1) )

!    Coastal Tiling
     FRAC_LAND  => D1(JFRAC_LAND:JFRAC_LAND+field_length(land_points,no_halo,1))
     TSTAR_LAND => D1(JTSTAR_LAND:JTSTAR_LAND+ &
      field_length(theta_points,no_halo,1))
     TSTAR_SEA  => D1(JTSTAR_SEA:JTSTAR_SEA+ &
      field_length(theta_points,no_halo,1))
     TSTAR_SICE => D1(JTSTAR_SICE:JTSTAR_SICE+ &
      field_length(theta_points,no_halo,1))

!    SeaIce and Land albedos
     SICE_ALB => D1(JSICE_ALB : JSICE_ALB+field_length(theta_points,no_halo,1))
     LAND_ALB => D1(JLAND_ALB : JLAND_ALB+field_length(theta_points,no_halo,1))

!    Large-Scale hydrology
     TI_MEAN   => D1(JTI_MEAN:JTI_MEAN+field_length(land_points,no_halo,1))
     TI_SIG    => D1(JTI_SIG:JTI_SIG+field_length(land_points,no_halo,1))
     FEXP      => D1(JFEXP:JFEXP+field_length(land_points,no_halo,1))
     GAMMA_INT => D1(JGAMMA_INT:JGAMMA_INT+field_length(land_points,no_halo,1))
     FSFC_SAT  => D1(JFSFC_SAT:JFSFC_SAT+field_length(land_points,no_halo,1))
     F_WETLAND => D1(JF_WETLAND:JF_WETLAND+field_length(land_points,no_halo,1))
     WATER_TABLE => D1(JWATER_TABLE:JWATER_TABLE+ &
      field_length(land_points,no_halo,1))

     STHZW   => D1(JSTHZW  : JSTHZW  + field_length(land_points,no_halo,1) )
     A_FSAT  => D1(JA_FSAT : JA_FSAT + field_length(land_points,no_halo,1) )
     C_FSAT  => D1(JC_FSAT : JC_FSAT + field_length(land_points,no_halo,1) )
     A_FWET  => D1(JA_FWET : JA_FWET + field_length(land_points,no_halo,1) )
     C_FWET  => D1(JC_FWET : JC_FWET + field_length(land_points,no_halo,1) )

!    Optional atmospheric primary variables
     ZH        => D1(JZH : JZH+field_length(theta_points,no_halo,1) )
     U_ADV     => D1(JU_ADV(1) : JU_ADV(1)+ &
      field_length(u_points,extended_halo,model_levels) )
     V_ADV     => D1(JV_ADV(1) : JV_ADV(1)+ &
      field_length(v_points,extended_halo,model_levels) )
     W_ADV     => D1(JW_ADV(0) : JW_ADV(0)+ &
      field_length(theta_points,extended_halo,model_levels+1) )
     NTML      => D1(JNTML : JNTML+field_length(theta_points,no_halo,1) )
     NBDSC     => D1(JNBDSC : JNBDSC+field_length(theta_points,no_halo,1) )
     NTDSC     => D1(JNTDSC : JNTDSC+field_length(theta_points,no_halo,1) )
     CUMULUS   => D1(JCUMULUS : JCUMULUS+field_length(theta_points,no_halo,1) )
     T1_SD     => D1(JT1_SD : JT1_SD+field_length(theta_points,no_halo,1) )
     Q1_SD     => D1(JQ1_SD : JQ1_SD+field_length(theta_points,no_halo,1) )
     CF_AREA   => D1(JCF_AREA(1) : JCF_AREA(1)+ &
      field_length(theta_points,no_halo,wet_levels) )
     CF_BULK   => D1(JCF_BULK(1) : JCF_BULK(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     CF_LIQUID => D1(JCF_LIQUID(1) : JCF_LIQUID(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     CF_FROZEN => D1(JCF_FROZEN(1) : JCF_FROZEN(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     ! size of cca varies according to L_3D_CCA
     CCA => D1(JCCA(1):JCCA(1)+field_length(theta_points,no_halo,n_cca_lev))
     CCB          => D1(JCCB : JCCB+field_length(theta_points,no_halo,1) )
     CCT          => D1(JCCT : JCCT+field_length(theta_points,no_halo,1) )
     CCLWP        => D1(JCCLWP : JCCLWP+field_length(theta_points,no_halo,1) )
     CANOPY_WATER => D1(JCANOPY_WATER : JCANOPY_WATER+ &
      field_length(land_points,no_halo,1) )
     LCBASE       => D1(JLCBASE : JLCBASE+field_length(theta_points,no_halo,1))
     CCW_RAD      => D1(JCCW_RAD(1) : JCCW_RAD(wet_levels)+theta_field_size)

!    Secondary Fields in D1
     EXNER_THETA_LEVELS => D1(JEXNER_THETA_LEVELS(1):JEXNER_THETA_LEVELS(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     P => D1(JP(1):JP(1)+field_length(theta_points,single_halo,model_levels+1))
     P_THETA_LEVELS => D1(JP_THETA_LEVELS(1):JP_THETA_LEVELS(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     PSTAR => D1(JPSTAR : JPSTAR+field_length(theta_points,no_halo,1) )
     SW_INCS => D1(JSW_INCS(0) : JSW_INCS(MODEL_LEVELS+1)+theta_field_size)
     LW_INCS => D1(JLW_INCS(0) : JLW_INCS(MODEL_LEVELS)+theta_field_size)

!    Direct PAR flux for STOCHEM
     DIRPAR => D1(JDIRPAR : JDIRPAR+field_length(theta_points,no_halo,1) )

!    Soil Fields
     SMCL => D1(JSMCL(1):JSMCL(1)+field_length(land_points,no_halo,sm_levels))
     DEEP_SOIL_TEMP => D1(J_DEEP_SOIL_TEMP(1):J_DEEP_SOIL_TEMP(1)+ &
      field_length(land_points,no_halo,sm_levels))
     VOL_SMC_WILT   => D1(JVOL_SMC_WILT:JVOL_SMC_WILT+ &
      field_length(land_points,no_halo,1))
     VOL_SMC_CRIT   => D1(JVOL_SMC_CRIT:JVOL_SMC_CRIT+ &
      field_length(land_points,no_halo,1))
     VOL_SMC_SAT    => D1(JVOL_SMC_SAT:JVOL_SMC_SAT+ &
      field_length(land_points,no_halo,1))
     SAT_SOIL_COND  => D1(JSAT_SOIL_COND:JSAT_SOIL_COND+ &
      field_length(land_points,no_halo,1))
     THERM_CAP  =>D1(JTHERM_CAP:JTHERM_CAP+field_length(land_points,no_halo,1))
     THERM_COND =>D1(JTHERM_COND:JTHERM_COND+ &
      field_length(land_points,no_halo,1))
     CLAPP_HORN =>D1(JCLAPP_HORN:JCLAPP_HORN+ &
      field_length(land_points,no_halo,1))
     SAT_SOILW_SUCTION => D1(JSAT_SOILW_SUCTION:JSAT_SOILW_SUCTION+ &
      field_length(land_points,no_halo,1))
     STHU => D1(JSTHU(1):JSTHU(1)+field_length(land_points,no_halo,sm_levels))
     STHF => D1(JSTHF(1):JSTHF(1)+field_length(land_points,no_halo,sm_levels))

!    Vegetation Fields
     Z0          => D1(JZ0:JZ0+field_length(theta_points,no_halo,1))
     VEG_FRAC    => D1(JVEG_FRAC:JVEG_FRAC+field_length(land_points,no_halo,1))
     ROOT_DEPTH  => D1(JROOT_DEPTH:JROOT_DEPTH+ &
      field_length(land_points,no_halo,1))
     SFA         => D1(JSFA:JSFA+field_length(land_points,no_halo,1))
     MDSA        => D1(JMDSA:JMDSA+field_length(land_points,no_halo,1))
     SURF_RESIST => D1(JSURF_RESIST:JSURF_RESIST+ &
      field_length(land_points,no_halo,1))
     SURF_CAP    => D1(JSURF_CAP:JSURF_CAP+field_length(land_points,no_halo,1))
     INFILT      => D1(JINFILT:JINFILT+field_length(land_points,no_halo,1))
     LAI         => D1(JLAI:JLAI+field_length(land_points,no_halo,1))
     CANHT       => D1(JCANHT:JCANHT+field_length(land_points,no_halo,1))
     GS          => D1(JGS:JGS+field_length(land_points,no_halo,1))

!    CABLE
     TSOIL_TILE  => D1(JTSOIL_TILE(1):JTSOIL_TILE(1)+  &
                         field_length(land_points,no_halo,sm_levels*ntiles))
     SMCL_TILE   => D1(JSMCL_TILE(1):JSMCL_TILE(1)+  &
                          field_length(land_points,no_halo,sm_levels*ntiles))
!  Not needed MRD
!!$     STHU_TILE   => D1(JSTHU_TILE:JSTHU_TILE+  &
!!$                          field_length(land_points,no_halo,sm_levels*ntiles))
     STHF_TILE   => D1(JSTHF_TILE(1):JSTHF_TILE(1)+  &
                   field_length(land_points,no_halo,sm_levels*ntiles))
! MRD - should be a parameter for number of snow levels here rather than 3
     SNOW_DEPTH3L=> D1(JSNOW_DEPTH3L(1):JSNOW_DEPTH3L(1)+  &
                          field_length(land_points,no_halo,3*ntiles))
     SNOW_MASS3L => D1(JSNOW_MASS3L(1):JSNOW_MASS3L(1)+  &
                          field_length(land_points,no_halo,3*ntiles))
     SNOW_TMP3L  => D1(JSNOW_TMP3L(1):JSNOW_TMP3L(1)+  &
                          field_length(land_points,no_halo,3*ntiles))
     SNOW_RHO3L  => D1(JSNOW_RHO3L(1):JSNOW_RHO3L(1)+  &
                          field_length(land_points,no_halo,3*ntiles))
     SNOW_RHO1L  => D1(JSNOW_RHO1L:JSNOW_RHO1L+  &
                          field_length(land_points,no_halo,ntiles))
     SNOW_AGE    => D1(JSNOW_AGE:JSNOW_AGE+  &
                          field_length(land_points,no_halo,ntiles))
     SNOW_FLG3L  => D1(JSNOW_FLG3L:JSNOW_FLG3L+  &
                          field_length(land_points,no_halo,ntiles))
! Lestevens Sept 2012 - CASA-CNP
     CPOOL_TILE  => D1(JCPOOL_TILE(1):JCPOOL_TILE(1)+ &
                          field_length(land_points,no_halo,10*ntiles))
     NPOOL_TILE  => D1(JNPOOL_TILE(1):JNPOOL_TILE(1)+ &
                          field_length(land_points,no_halo,10*ntiles))
     PPOOL_TILE  => D1(JPPOOL_TILE(1):JPPOOL_TILE(1)+ &
                          field_length(land_points,no_halo,12*ntiles))
     SOIL_ORDER  => D1(JSOIL_ORDER:JSOIL_ORDER+ &
                          field_length(land_points,no_halo,1))
     NIDEP       => D1(JNIDEP:JNIDEP+ &
                          field_length(land_points,no_halo,1))
     NIFIX       => D1(JNIFIX:JNIFIX+ &
                          field_length(land_points,no_halo,1))
     PWEA        => D1(JPWEA:JPWEA+ &
                          field_length(land_points,no_halo,1))
     PDUST       => D1(JPDUST:JPDUST+ &
                          field_length(land_points,no_halo,1))
     GLAI        => D1(JGLAI:JGLAI+ &
                          field_length(land_points,no_halo,ntiles))
     PHENPHASE   => D1(JPHENPHASE:JPHENPHASE+ &
                          field_length(land_points,no_halo,ntiles))
     CABLE_LAI => D1(JCABLE_LAI(1):JCABLE_LAI(1)+  &
                         field_length(land_points,no_halo,ntiles)) 

!    Orography Fields
     OROG_SIL => D1(JOROG_SIL : JOROG_SIL+field_length(land_points,no_halo,1) )
     OROG_HO2 => D1(JOROG_HO2 : JOROG_HO2+field_length(land_points,no_halo,1) )
     OROG_SD => D1(JOROG_SD : JOROG_SD+field_length(land_points,no_halo,1) )
     OROG_GRAD_X  => D1(JOROG_GRAD_X : JOROG_GRAD_X+ &
      field_length(land_points,no_halo,1) )
     OROG_GRAD_Y  => D1(JOROG_GRAD_Y : JOROG_GRAD_Y+ &
      field_length(land_points,no_halo,1) )
     OROG_GRAD_XX => D1(JOROG_GRAD_XX : JOROG_GRAD_XX+ &
      field_length(land_points,no_halo,1) )
     OROG_GRAD_XY => D1(JOROG_GRAD_XY : JOROG_GRAD_XY+ &
      field_length(land_points,no_halo,1) )
     OROG_GRAD_YY => D1(JOROG_GRAD_YY : JOROG_GRAD_YY+ &
      field_length(land_points,no_halo,1) )

!    Sea/Sea Ice Fields
     U_SEA => D1(JU_SEA : JU_SEA+field_length(u_points,no_halo,1) )
     V_SEA => D1(JV_SEA : JV_SEA+field_length(v_points,no_halo,1) )
     ICE_FRACTION  => D1(JICE_FRACTION  : JICE_FRACTION + &
      field_length(theta_points_sea_only,no_halo,1) )
     ICE_THICKNESS => D1(JICE_THICKNESS : JICE_THICKNESS+ &
      field_length(theta_points_sea_only,no_halo,1) )
     TI => D1(JTI : JTI+field_length(theta_points_sea_only,no_halo,1) )
     ICE_FRACT_CAT => D1(JICE_FRACT_CAT : JICE_FRACT_CAT+ &
      field_length(theta_points,no_halo,nice) )
     ICE_THICK_CAT => D1(JICE_THICK_CAT : JICE_THICK_CAT+ &
      field_length(theta_points,no_halo,nice) )
     TI_CAT => D1(JTI_CAT : JTI_CAT+field_length(theta_points,no_halo,nice) )
     U_0_P => D1(JU_0_P : JU_0_P+field_length(theta_points,no_halo,1) )
     V_0_P => D1(JV_0_P : JV_0_P+field_length(theta_points,no_halo,1) )

!    Snow Fields
     SNODEP => D1(JSNODEP : JSNODEP+field_length(theta_points,no_halo,1) )
     SNODEP_SEA => D1(JSNODEP_SEA : JSNODEP_SEA+ &
      field_length(theta_points_sea_only,no_halo,1) )
     SNODEP_SEA_CAT => D1(JSNODEP_SEA_CAT : JSNODEP_SEA_CAT+ &
      field_length(theta_points,no_halo,nice) )
! SNSOOT may not be used as of vn6.6
     SNSOOT => D1(JSNSOOT : JSNSOOT+field_length(theta_points,no_halo,1) )
     CATCH_SNOW => D1(JCATCH_SNOW : JCATCH_SNOW+ &
      field_length(land_points,no_halo,1) )
     SNOW_GRND => D1(JSNOW_GRND : JSNOW_GRND+ &
      field_length(land_points,no_halo,1) )

!    OZONE
     O3 => D1(JOZONE(1) : JOZONE(1)+ &
      field_length(ozone_points,no_halo,ozone_levels) )
!    Tropopause-based Ozone
     IF (TPPS_OZONE_LEVELS > 0) THEN
       TPPSOZONE => D1(JTPPSOZONE(1) : JTPPSOZONE(1)+ &
        field_length(ozone_points,no_halo,ozone_levels) )
     ELSE
       TPPSOZONE => dummy_field
     END IF

!    Ozone tracer field and cariolle parameters
     OZONE_TRACER => D1(JOZONE_TRACER(1)   : JOZONE_TRACER(1)+       &
                 field_length(theta_points,single_halo,model_levels))
     O3_PROD_LOSS => D1(JO3_PROD_LOSS(1)   : JO3_PROD_LOSS(1)+       &
                             (ROWS*MODEL_LEVELS))
     O3_P_L_VMR   => D1(JO3_P_L_VMR(1)   : JO3_P_L_VMR(1)+           &
                             (ROWS*MODEL_LEVELS))
     O3_VMR       => D1(JO3_VMR (1)   : JO3_VMR(1)+                  &
                             (ROWS*MODEL_LEVELS))
     O3_P_L_TEMP  => D1(JO3_P_L_TEMP(1)   : JO3_P_L_TEMP(1)+         &
                             (ROWS*MODEL_LEVELS))
     O3_TEMP      => D1(JO3_TEMP(1)   : JO3_TEMP(1)+                 &
                             (ROWS*MODEL_LEVELS))    
     O3_P_L_COLO3 => D1(JO3_P_L_COLO3(1)   : JO3_P_L_COLO3(1)+       &
                             (ROWS*MODEL_LEVELS))
     O3_COLO3     => D1(JO3_COLO3(1)   : JO3_COLO3(1)+               &
                             (ROWS*MODEL_LEVELS))
!    STOCHEM fields
     CH4_STOCH => D1(JCH4_STOCH(1) : JCH4_STOCH(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     O3_STOCH  => D1(JO3_STOCH(1)  : JO3_STOCH(1) + &
      field_length(theta_points,no_halo,model_levels) )

!    Sources and Aerosol Ancillaries
     MURK_SOURCE => D1(JMURK_SOURCE(1) : JMURK_SOURCE(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     SO2_EM      => D1(JSO2_EM : JSO2_EM+field_length(theta_points,no_halo,1) )
     DMS_EM      => D1(JDMS_EM : JDMS_EM+field_length(theta_points,no_halo,1) )
     MURK        => D1(JMURK(1) : JMURK(1)+ &
      field_length(theta_points,single_halo,model_levels) )

!    Sulphur cycle
     SO2         => D1(JSO2(1) : JSO2(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DMS         => D1(JDMS(1) : JDMS(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SO4_AITKEN  => D1(JSO4_AITKEN(1) : JSO4_AITKEN(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SO4_ACCU    => D1(JSO4_ACCU(1) : JSO4_ACCU(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SO4_DISS    => D1(JSO4_DISS(1) : JSO4_DISS(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     H2O2        => D1(JH2O2(1) : JH2O2(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     NH3         => D1(JNH3(1) : JNH3(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SOOT_NEW    => D1(JSOOT_NEW(1) : JSOOT_NEW(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SOOT_AGD    => D1(JSOOT_AGD(1) : JSOOT_AGD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SOOT_CLD    => D1(JSOOT_CLD(1) : JSOOT_CLD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     BMASS_NEW   => D1(JBMASS_NEW(1) : JBMASS_NEW(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     BMASS_AGD   => D1(JBMASS_AGD(1) : JBMASS_AGD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     BMASS_CLD   => D1(JBMASS_CLD(1) : JBMASS_CLD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     OCFF_NEW    => D1(JOCFF_NEW(1) : JOCFF_NEW(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     OCFF_AGD    => D1(JOCFF_AGD(1) : JOCFF_AGD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     OCFF_CLD    => D1(JOCFF_CLD(1) : JOCFF_CLD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SO2_NATEM   => D1(JSO2_NATEM(1) : JSO2_NATEM(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     OH          => D1(JOH(1) : JOH(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     HO2         => D1(JHO2(1) : JHO2(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     H2O2_LIMIT  => D1(JH2O2_LIMIT(1) : JH2O2_LIMIT(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     O3_CHEM     => D1(JO3_CHEM(1) : JO3_CHEM(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     SO2_HILEM   => D1(JSO2_HILEM : JSO2_HILEM+ &
      field_length(theta_points,no_halo,1) )
     NH3_EM      => D1(JNH3_EM : JNH3_EM+ &
      field_length(theta_points,no_halo,1) )
     SOOT_EM     => D1(JSOOT_EM : JSOOT_EM+ &
      field_length(theta_points,no_halo,1) )
     SOOT_HILEM  => D1(JSOOT_HILEM : JSOOT_HILEM+ &
      field_length(theta_points,no_halo,1) )
     BMASS_EM    => D1(JBMASS_EM : JBMASS_EM+ &
      field_length(theta_points,no_halo,1) )
     BMASS_HILEM => D1(JBMASS_HILEM : JBMASS_HILEM+ &
      field_length(theta_points,no_halo,1) )
     OCFF_EM     => D1(JOCFF_EM : JOCFF_EM+ &
      field_length(theta_points,no_halo,1) )
     OCFF_HILEM  => D1(JOCFF_HILEM : JOCFF_HILEM+ &
      field_length(theta_points,no_halo,1) )
     DMS_CONC    => D1(JDMS_CONC : JDMS_CONC+ &
      field_length(theta_points,no_halo,1) )
     DMS_OFLUX   => D1(JDMS_OFLUX : JDMS_OFLUX+ &
      field_length(theta_points,no_halo,1) )

! Aerosol climatologies
     ARCLBIOG_BG => D1(JARCLBIOG_BG(1) : JARCLBIOG_BG(model_levels)+theta_field_size)
     ARCLBIOM_FR => D1(JARCLBIOM_FR(1) : JARCLBIOM_FR(model_levels)+theta_field_size)
     ARCLBIOM_AG => D1(JARCLBIOM_AG(1) : JARCLBIOM_AG(model_levels)+theta_field_size)
     ARCLBIOM_IC => D1(JARCLBIOM_IC(1) : JARCLBIOM_IC(model_levels)+theta_field_size)
     ARCLBLCK_FR => D1(JARCLBLCK_FR(1) : JARCLBLCK_FR(model_levels)+theta_field_size)
     ARCLBLCK_AG => D1(JARCLBLCK_AG(1) : JARCLBLCK_AG(model_levels)+theta_field_size)
     ARCLSSLT_FI => D1(JARCLSSLT_FI(1) : JARCLSSLT_FI(model_levels)+theta_field_size)
     ARCLSSLT_JT => D1(JARCLSSLT_JT(1) : JARCLSSLT_JT(model_levels)+theta_field_size)
     ARCLSULP_AC => D1(JARCLSULP_AC(1) : JARCLSULP_AC(model_levels)+theta_field_size)
     ARCLSULP_AK => D1(JARCLSULP_AK(1) : JARCLSULP_AK(model_levels)+theta_field_size)
     ARCLSULP_DI => D1(JARCLSULP_DI(1) : JARCLSULP_DI(model_levels)+theta_field_size)
     ARCLDUST_B1 => D1(JARCLDUST_B1(1) : JARCLDUST_B1(model_levels)+theta_field_size)
     ARCLDUST_B2 => D1(JARCLDUST_B2(1) : JARCLDUST_B2(model_levels)+theta_field_size)
     ARCLDUST_B3 => D1(JARCLDUST_B3(1) : JARCLDUST_B3(model_levels)+theta_field_size)
     ARCLDUST_B4 => D1(JARCLDUST_B4(1) : JARCLDUST_B4(model_levels)+theta_field_size)
     ARCLDUST_B5 => D1(JARCLDUST_B5(1) : JARCLDUST_B5(model_levels)+theta_field_size)
     ARCLDUST_B6 => D1(JARCLDUST_B6(1) : JARCLDUST_B6(model_levels)+theta_field_size)
     ARCLOCFF_FR => D1(JARCLOCFF_FR(1) : JARCLOCFF_FR(model_levels)+theta_field_size)
     ARCLOCFF_AG => D1(JARCLOCFF_AG(1) : JARCLOCFF_AG(model_levels)+theta_field_size)
     ARCLOCFF_IC => D1(JARCLOCFF_IC(1) : JARCLOCFF_IC(model_levels)+theta_field_size)
     ARCLDLTA_DL => D1(JARCLDLTA_DL(1) : JARCLDLTA_DL(model_levels)+theta_field_size)
     
!    Mineral Dust Schema
     SOIL_CLAY  => D1(JSOIL_CLAY : JSOIL_CLAY+ &
      field_length(theta_points,no_halo,1) )
     SOIL_SILT  => D1(JSOIL_SILT : JSOIL_SILT+ &
      field_length(theta_points,no_halo,1) )
     SOIL_SAND  => D1(JSOIL_SAND : JSOIL_SAND+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL1 => D1(JDUST_MREL1 : JDUST_MREL1+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL2 => D1(JDUST_MREL2 : JDUST_MREL2+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL3 => D1(JDUST_MREL3 : JDUST_MREL3+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL4 => D1(JDUST_MREL4 : JDUST_MREL4+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL5 => D1(JDUST_MREL5 : JDUST_MREL5+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL6 => D1(JDUST_MREL6 : JDUST_MREL6+ &
      field_length(theta_points,no_halo,1) )
     DUST_DIV1  => D1(JDUST_DIV1(1) : JDUST_DIV1(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DUST_DIV2  => D1(JDUST_DIV2(1) : JDUST_DIV2(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DUST_DIV3  => D1(JDUST_DIV3(1) : JDUST_DIV3(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DUST_DIV4  => D1(JDUST_DIV4(1) : JDUST_DIV4(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DUST_DIV5  => D1(JDUST_DIV5(1) : JDUST_DIV5(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DUST_DIV6  => D1(JDUST_DIV6(1) : JDUST_DIV6(1)+ &
      field_length(theta_points,single_halo,model_levels) )

!    Carbon Cycle
     CO2FLUX   => D1(J_CO2FLUX : J_CO2FLUX+ &
      field_length(theta_points,no_halo,1) )
     CO2_EMITS => D1(J_CO2_EMITS : J_CO2_EMITS+ &
      field_length(theta_points,no_halo,1) )
     CO2       => D1(JCO2(1) : JCO2(model_levels)+theta_off_size)

!    level dependent constants
     zseak_theta => D1(JZSEAK_THETA : JZSEAK_THETA+(model_levels+1) )
     Ck_theta    => D1(JCK_THETA : JCK_THETA+(model_levels+1) )
     zseak_rho   => D1(JZSEAK_RHO : JZSEAK_RHO+(model_levels+1) )
     Ck_rho      => D1(JCK_RHO : JCK_RHO+(model_levels+1) )

!    Tracer Fluxes - kdcorbin, 05/10
     TRACER_FLUX1  => D1(JTRACER_FLUX1 : JTRACER_FLUX1 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX2  => D1(JTRACER_FLUX2 : JTRACER_FLUX2 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX3  => D1(JTRACER_FLUX3 : JTRACER_FLUX3 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX4  => D1(JTRACER_FLUX4 : JTRACER_FLUX4 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX5  => D1(JTRACER_FLUX5 : JTRACER_FLUX5 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX6  => D1(JTRACER_FLUX6 : JTRACER_FLUX6 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX7  => D1(JTRACER_FLUX7 : JTRACER_FLUX7 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX8  => D1(JTRACER_FLUX8 : JTRACER_FLUX8 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX9  => D1(JTRACER_FLUX9 : JTRACER_FLUX9 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX10 => D1(JTRACER_FLUX10: JTRACER_FLUX10 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX11 => D1(JTRACER_FLUX11: JTRACER_FLUX11 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX12 => D1(JTRACER_FLUX12: JTRACER_FLUX12 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX13 => D1(JTRACER_FLUX13: JTRACER_FLUX13 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX14 => D1(JTRACER_FLUX14: JTRACER_FLUX14 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX15 => D1(JTRACER_FLUX15: JTRACER_FLUX15 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX16 => D1(JTRACER_FLUX16: JTRACER_FLUX16 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX17 => D1(JTRACER_FLUX17: JTRACER_FLUX17 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX18 => D1(JTRACER_FLUX18: JTRACER_FLUX18 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX19 => D1(JTRACER_FLUX19: JTRACER_FLUX19 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX20 => D1(JTRACER_FLUX20: JTRACER_FLUX20 + &
      field_length(theta_points,no_halo,1) )

!    User ancillaries
     USER_ANC1   => D1(JUSER_ANC1  : JUSER_ANC1 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC2   => D1(JUSER_ANC2  : JUSER_ANC2 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC3   => D1(JUSER_ANC3  : JUSER_ANC3 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC4   => D1(JUSER_ANC4  : JUSER_ANC4 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC5   => D1(JUSER_ANC5  : JUSER_ANC5 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC6   => D1(JUSER_ANC6  : JUSER_ANC6 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC7   => D1(JUSER_ANC7  : JUSER_ANC7 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC8   => D1(JUSER_ANC8  : JUSER_ANC8 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC9   => D1(JUSER_ANC9  : JUSER_ANC9 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC10  => D1(JUSER_ANC10 : JUSER_ANC10+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC11  => D1(JUSER_ANC11 : JUSER_ANC11+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC12  => D1(JUSER_ANC12 : JUSER_ANC12+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC13  => D1(JUSER_ANC13 : JUSER_ANC13+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC14  => D1(JUSER_ANC14 : JUSER_ANC14+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC15  => D1(JUSER_ANC15 : JUSER_ANC15+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC16  => D1(JUSER_ANC16 : JUSER_ANC16+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC17  => D1(JUSER_ANC17 : JUSER_ANC17+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC18  => D1(JUSER_ANC18 : JUSER_ANC18+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC19  => D1(JUSER_ANC19 : JUSER_ANC19+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC20  => D1(JUSER_ANC20 : JUSER_ANC20+ &
      field_length(theta_points,no_halo,1) )
     USER_MULT1  => D1(JUSER_MULT1(1)  : JUSER_MULT1(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT2  => D1(JUSER_MULT2(1)  : JUSER_MULT2(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT3  => D1(JUSER_MULT3(1)  : JUSER_MULT3(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT4  => D1(JUSER_MULT4(1)  : JUSER_MULT4(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT5  => D1(JUSER_MULT5(1)  : JUSER_MULT5(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT6  => D1(JUSER_MULT6(1)  : JUSER_MULT6(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT7  => D1(JUSER_MULT7(1)  : JUSER_MULT7(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT8  => D1(JUSER_MULT8(1)  : JUSER_MULT8(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT9  => D1(JUSER_MULT9(1)  : JUSER_MULT9(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT10 => D1(JUSER_MULT10(1) : JUSER_MULT10(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT11 => D1(JUSER_MULT11(1) : JUSER_MULT11(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT12 => D1(JUSER_MULT12(1) : JUSER_MULT12(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT13 => D1(JUSER_MULT13(1) : JUSER_MULT13(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT14 => D1(JUSER_MULT14(1) : JUSER_MULT14(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT15 => D1(JUSER_MULT15(1) : JUSER_MULT15(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT16 => D1(JUSER_MULT16(1) : JUSER_MULT16(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT17 => D1(JUSER_MULT17(1) : JUSER_MULT17(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT18 => D1(JUSER_MULT18(1) : JUSER_MULT18(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT19 => D1(JUSER_MULT19(1) : JUSER_MULT19(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT20 => D1(JUSER_MULT20(1) : JUSER_MULT20(1)+ &
      field_length(theta_points,no_halo,model_levels) )

!    Tiled vegetation and triffid
     FRAC_TYP  => D1(JFRAC_TYP:JFRAC_TYP+field_length(land_points,no_halo,1))
     FRAC_CON1 => D1(JFRAC_CON1:JFRAC_CON1+field_length(land_points,no_halo,1))
     FRAC_CON2 => D1(JFRAC_CON2:JFRAC_CON2+field_length(land_points,no_halo,1))
     FRAC_CON3 => D1(JFRAC_CON3:JFRAC_CON3+field_length(land_points,no_halo,1))
     FRAC_CON4 => D1(JFRAC_CON4:JFRAC_CON4+field_length(land_points,no_halo,1))
     FRAC_CON5 => D1(JFRAC_CON5:JFRAC_CON5+field_length(land_points,no_halo,1))
     FRAC_CON6 => D1(JFRAC_CON6:JFRAC_CON6+field_length(land_points,no_halo,1))
     FRAC_CON7 => D1(JFRAC_CON7:JFRAC_CON7+field_length(land_points,no_halo,1))
     FRAC_CON8 => D1(JFRAC_CON8:JFRAC_CON8+field_length(land_points,no_halo,1))
     FRAC_CON9 => D1(JFRAC_CON9:JFRAC_CON9+field_length(land_points,no_halo,1))
     LAI_PFT   => D1(JLAI_PFT:JLAI_PFT+field_length(land_points,no_halo,1))
     CANHT_PFT => D1(JCANHT_PFT:JCANHT_PFT+field_length(land_points,no_halo,1))
     DISTURB_VEG => D1(JDISTURB:JDISTURB+field_length(land_points,no_halo,1))
     SOIL_ALB  => D1(JSOIL_ALB:JSOIL_ALB+field_length(land_points,no_halo,1))
     SOIL_CARB => D1(JSOIL_CARB:JSOIL_CARB+field_length(land_points,no_halo,1))
     SOIL_CARB1 => D1(JSOIL_CARB1:JSOIL_CARB1+ &
      field_length(land_points,no_halo,1))
     SOIL_CARB2 => D1(JSOIL_CARB2:JSOIL_CARB2+ &
      field_length(land_points,no_halo,1))
     SOIL_CARB3 => D1(JSOIL_CARB3:JSOIL_CARB3+ &
      field_length(land_points,no_halo,1))
     SOIL_CARB4 => D1(JSOIL_CARB4:JSOIL_CARB4+ &
      field_length(land_points,no_halo,1))
     NPP_PFT_ACC    => D1(JNPP_PFT_ACC:JNPP_PFT_ACC+ &
      field_length(land_points,no_halo,1))
     G_LF_PFT_ACC   => D1(JG_LF_PFT_ACC:JG_LF_PFT_ACC+ &
      field_length(land_points,no_halo,1))
     G_PHLF_PFT_ACC => D1(JG_PHLF_PFT_ACC:JG_PHLF_PFT_ACC+ &
      field_length(land_points,no_halo,1))
     RSP_W_PFT_ACC  => D1(JRSP_W_PFT_ACC:JRSP_W_PFT_ACC+ &
      field_length(land_points,no_halo,1))
     RSP_S_ACC  => D1(JRSP_S_ACC:JRSP_S_ACC+ &
      field_length(land_points,no_halo,1))
     RSP_S_ACC1 => D1(JRSP_S_ACC1:JRSP_S_ACC1+ &
      field_length(land_points,no_halo,1))
     RSP_S_ACC2 => D1(JRSP_S_ACC2:JRSP_S_ACC2+ &
      field_length(land_points,no_halo,1))
     RSP_S_ACC3 => D1(JRSP_S_ACC3:JRSP_S_ACC3+ &
      field_length(land_points,no_halo,1))
     RSP_S_ACC4 => D1(JRSP_S_ACC4:JRSP_S_ACC4+ &
      field_length(land_points,no_halo,1))
     CAN_WATER_TILE => D1(JCAN_WATER_TILE:JCAN_WATER_TILE+ &
      field_length(land_points,no_halo,1))
     CATCH_TILE  => D1(JCATCH_TILE:JCATCH_TILE+ &
      field_length(land_points,no_halo,1))
     RGRAIN_TILE => D1(JRGRAIN_TILE:JRGRAIN_TILE+ &
      field_length(land_points,no_halo,1))
     ! TSNOW => no longer used
     TSTAR_TILE  => D1(JTSTAR_TILE:JTSTAR_TILE+ &
      field_length(land_points,no_halo,1))
     Z0_TILE     => D1(JZ0_TILE:JZ0_TILE+field_length(land_points,no_halo,1))
     SNODEP_TILE => D1(JSNODEP_TILE:JSNODEP_TILE+ &
      field_length(land_points,no_halo,1))
     INFIL_TILE  => D1(JINFIL_TILE:JINFIL_TILE+ &
      field_length(land_points,no_halo,1))
     DOLR_FIELD  => D1(JDOLR:JDOLR+field_length(theta_points,no_halo,1))
     LW_DOWN     => D1(JLW_DOWN:JLW_DOWN+field_length(theta_points,no_halo,1))
     SW_TILE_RTS => D1(JSW_TILE:JSW_TILE+field_length(land_points,no_halo,1))

!    River routing fields
     RIV_SEQUENCE  => D1(JRIV_SEQUENCE : JRIV_SEQUENCE+ &
      field_length(river_points,no_halo,1) )
     RIV_DIRECTION => D1(JRIV_DIRECTION : JRIV_DIRECTION+ &
      field_length(river_points,no_halo,1) )
     RIV_STORAGE   => D1(JRIV_STORAGE : JRIV_STORAGE+ &
      field_length(river_points,no_halo,1) )
     TOT_SURFROFF  => D1(JTOT_SURFROFF : JTOT_SURFROFF+ &
      field_length(river_points,no_halo,1) )
     TOT_SUBROFF   => D1(JTOT_SUBROFF : JTOT_SUBROFF+ &
      field_length(river_points,no_halo,1) )
     RIV_INLANDATM => D1(JRIV_INLANDATM : JRIV_INLANDATM+ &
      field_length(land_points,no_halo,1) )
     ! these are unitialised upon entering ATM_STEP
     RIV_IAREA     => dummy_field !D1(1:1+row_length*rows)
     RIV_SLOPE     => dummy_field !D1(1:1+row_length*rows)
     RIV_FLOWOBS1  => dummy_field !D1(1:1+row_length*rows)
     RIV_INEXT     => dummy_field !D1(1:1+row_length*rows)
     RIV_JNEXT     => dummy_field !D1(1:1+row_length*rows)
     RIV_LAND      => dummy_field !D1(1:1+row_length*rows)
     RIV_SUBSTORE  => dummy_field !D1(1:1+row_length*rows)
     RIV_SURFSTORE => dummy_field !D1(1:1+row_length*rows)
     RIV_FLOWIN    => dummy_field !D1(1:1+row_length*rows)
     RIV_BFLOWIN   => dummy_field !D1(1:1+row_length*rows)

!    Fields to be tretained in dumps for coupled models using OASIS

     C_DOWN =>  D1(JC_DOWN : JC_DOWN + &
                    field_length(theta_points,no_halo,1))

 
! gol124: auscom coupling
! probably not needed: coupling code seems to use def
! in include/argument/arg_atm_fields.h

!    Required for energy correction
     NET_FLUX  => D1(JNET_FLUX:JNET_FLUX+field_length(theta_points,no_halo,1))
     NET_MFLUX => D1(JNET_MFLUX:JNET_MFLUX+field_length(theta_points,no_halo,1))

!    Fields carried forward from previous version
     TSTAR_ANOM => D1(JTSTAR_ANOM : JTSTAR_ANOM+field_length(theta_points,no_halo,1) )

! SLAB model removed as of vn7.0
!#if defined (SLAB)
!
!     TCLIM => D1(JTCLIM : JTCLIM+field_length(theta_points,no_halo,1))
!     HCLIM => D1(JHCLIM : JHCLIM+field_length(theta_points,no_halo,1))
!     TSLAB => D1(JTSLAB : JTSLAB+field_length(theta_points_sea_only,no_halo,1))
!     CHEAT => D1(JCHEAT : JCHEAT+field_length(theta_points_sea_only,no_halo,1))
!     OIFLX => D1(JOIFLX : JOIFLX+field_length(theta_points_sea_only,no_halo,1))
!     UICE  => D1(JUICE  : JUICE+field_length(u_points,no_halo,1))
!     VICE  => D1(JVICE  : JVICE+field_length(v_points,no_halo,1))
!     SIG11NE => D1(JSIG11NE : JSIG11NE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG11SE => D1(JSIG11SE : JSIG11SE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG11SW => D1(JSIG11SW : JSIG11SW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG11NW => D1(JSIG11NW : JSIG11NW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG12NE => D1(JSIG12NE : JSIG12NE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG12SE => D1(JSIG12SE : JSIG12SE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG12SW => D1(JSIG12SW : JSIG12SW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG12NW => D1(JSIG12NW : JSIG12NW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG22NE => D1(JSIG22NE : JSIG22NE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG22SE => D1(JSIG22SE : JSIG22SE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG22SW => D1(JSIG22SW : JSIG22SW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG22NW => D1(JSIG22NW : JSIG22NW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!
!#endif ! defined (SLAB)

!    lateral boundary conditions

     OROG_LBC  => D1(JOROG_LBC : JOROG_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*1) )
     U_LBC     => D1(JU_LBC : JU_LBC + &
      (LENRIMA(fld_type_u,halo_type_extended,1)*model_levels) )
     V_LBC     => D1(JV_LBC : JU_LBC +  &
     (LENRIMA(fld_type_v,halo_type_extended,1)*model_levels) )
     W_LBC     => D1(JW_LBC : JW_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,rima_type_norm)*(model_levels+1)) )
     RHO_LBC   => D1(JRHO_LBC : JRHO_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*model_levels) )
     THETA_LBC => D1(JTHETA_LBC : JTHETA_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*model_levels) )
     Q_LBC     => D1(JQ_LBC : JQ_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QCL_LBC   => D1(JQCL_LBC : JQCL_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QCF_LBC   => D1(JQCF_LBC : JQCF_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     EXNER_LBC => D1(JEXNER_LBC : JEXNER_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*(model_levels+1)) )
     U_ADV_LBC => D1(JU_ADV_LBC : JU_ADV_LBC + &
      (LENRIMA(fld_type_u,halo_type_extended,1)*model_levels) )
     V_ADV_LBC => D1(JV_ADV_LBC : JV_ADV_LBC + &
      (LENRIMA(fld_type_v,halo_type_extended,1)*model_levels) )
     W_ADV_LBC => D1(JW_ADV_LBC : JW_ADV_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*(model_levels+1)) )
     QCF2_LBC  => D1(JQCF2_LBC : JQCF2_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QRAIN_LBC => D1(JQRAIN_LBC : JQRAIN_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QGRAUP_LBC  => D1(JQGRAUP_LBC : JQGRAUP_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_BULK_LBC => D1(JCF_BULK_LBC : JCF_BULK_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_LIQUID_LBC => D1(JCF_LIQUID_LBC : JCF_LIQUID_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_FROZEN_LBC => D1(JCF_FROZEN_LBC : JCF_FROZEN_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     MURK_LBC  => D1(JMURK_LBC : JMURK_LBC + &
      (LENRIMA(fld_type_p,halo_type_single,1)*model_levels) )

     U_LBC_TEND     => D1(JU_LBC_TEND : JU_LBC_TEND + &
      (LENRIMA(fld_type_u,halo_type_extended,1)*model_levels) )
     V_LBC_TEND     => D1(JV_LBC_TEND : JV_LBC_TEND + &
      (LENRIMA(fld_type_v,halo_type_extended,1)*model_levels) )
     W_LBC_TEND     => D1(JW_LBC_TEND : JW_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,rima_type_norm)*(model_levels+1)) )
     RHO_LBC_TEND   => D1(JRHO_LBC_TEND : JRHO_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*model_levels) )
     THETA_LBC_TEND => D1(JTHETA_LBC_TEND : JTHETA_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*model_levels) )
     Q_LBC_TEND     => D1(JQ_LBC_TEND : JQ_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QCL_LBC_TEND   => D1(JQCL_LBC_TEND : JQCL_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QCF_LBC_TEND   => D1(JQCF_LBC_TEND : JQCF_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     EXNER_LBC_TEND => D1(JEXNER_LBC_TEND : JEXNER_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*(model_levels+1)) )
     U_ADV_LBC_TEND => D1(JU_ADV_LBC_TEND : JU_ADV_LBC_TEND + &
      (LENRIMA(fld_type_u,halo_type_extended,1)*model_levels) )
     V_ADV_LBC_TEND => D1(JV_ADV_LBC_TEND : JV_ADV_LBC_TEND + &
      (LENRIMA(fld_type_v,halo_type_extended,1)*model_levels) )
     W_ADV_LBC_TEND => D1(JW_ADV_LBC_TEND : JW_ADV_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*(model_levels+1)) )
     QCF2_LBC_TEND => D1(JQCF2_LBC_TEND : JQCF2_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QRAIN_LBC_TEND => D1(JQRAIN_LBC_TEND : JQRAIN_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QGRAUP_LBC_TEND => D1(JQGRAUP_LBC_TEND : JQGRAUP_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_BULK_LBC_TEND => D1(JCF_BULK_LBC_TEND : JCF_BULK_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_LIQUID_LBC_TEND => D1(JCF_LIQUID_LBC_TEND : JCF_LIQUID_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_FROZEN_LBC_TEND => D1(JCF_FROZEN_LBC_TEND : JCF_FROZEN_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     MURK_LBC_TEND => D1(JMURK_LBC_TEND : JMURK_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_single,1)*model_levels) )

! Oxidant concentrations from UKCA for use in HadGEM sulphur
! cycle (these are in Section 33):
      IF(L_SULPC_ONLINE_OXIDANTS .AND. L_UKCA) THEN 
         OH_UKCA   => D1(JOH_UKCA(1)   : JOH_UKCA(1)+ &
          field_length(theta_points,single_halo,model_levels))
         H2O2_UKCA => D1(JH2O2_UKCA(1) : JH2O2_UKCA(1)+ &
          field_length(theta_points,single_halo,model_levels))
         HO2_UKCA  => D1(JHO2_UKCA(1)  : JHO2_UKCA(1)+ &
          field_length(theta_points,single_halo,model_levels))
         O3_UKCA   => D1(JO3_UKCA(1)   : JO3_UKCA(1)+ &
          field_length(theta_points,single_halo,model_levels))
      ELSE 
        OH_UKCA   => dummy_field
        H2O2_UKCA => dummy_field
        HO2_UKCA  => dummy_field
        O3_UKCA   => dummy_field
      END IF  
! 1.1 point tracer fields to D1

      ! find out how many tracers are active
      nActiveTracers=0
      DO nTracer=A_TRACER_FIRST,A_TRACER_LAST
        IF (SI(nTracer,33,atmos_im) /= 1) THEN
          nActiveTracers = nActiveTracers+1
        END IF
      END DO ! nTracer
      IF (nActiveTracers /= 0) THEN
        ! set the pointer to the appropriate section of D1
        TRACER => D1(JTRACER(1,A_TRACER_FIRST) : &
         JTRACER(tr_levels,nActiveTracers)+theta_off_size)
        TRACER_LBC => D1(JTRACER_LBC(A_TRACER_FIRST) : &
         JTRACER_LBC(nActiveTracers) + &
         (LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels) )
        TRACER_LBC_TEND => D1(JTRACER_LBC_TEND(A_TRACER_FIRST) : &
         JTRACER_LBC_TEND(nActiveTracers) + &
         (LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels) )
      ELSE
        ! or set it to something non-null if there are no active tracers
        TRACER => dummy_field
        TRACER_LBC => dummy_field
        TRACER_LBC_TEND => dummy_field
      END IF
     
      ! do the same for section 34 (UKCA) tracers     
      nActiveTracers=0
      DO nTracer=A_UKCA_FIRST,A_UKCA_LAST
        IF (SI(nTracer,34,atmos_im) /= 1) THEN
          nActiveTracers = nActiveTracers+1
        END IF
      END DO ! nTracer          
      IF (nActiveTracers /= 0) THEN
        TRACER_UKCA => D1(JTR_UKCA(1,A_UKCA_FIRST) : &
         JTR_UKCA(tr_levels,nActiveTracers)+theta_off_size)
      ELSE
        TRACER_UKCA => dummy_field
      END IF

END SUBROUTINE Set_Atm_Fields

