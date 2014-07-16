
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine AC_CTL   -----------------------------------------------
!LL
!LL Level 2 control routine
!LL version for CRAY YMP
!LL
!LL programming standard : unified model documentation paper No 3
!LL
!LL Logical components covered : P1
!LL
!LL Project task : P0
!LLEND -----------------------------------------------------------------
!*L Arguments
      SUBROUTINE AC_CTL(INT18,P_FIELDDA,Q_LEVELSDA,                     &
     & model_levelsda, theta_star_inc_halo,                             &
     & q_star, qcl_star, qcf_star,                                      &
     & OBS_FLAG,OBS,obs_flag_len,obs_len,                               &
     & p, p_theta_levels, exner_theta_levels, r_theta_levels,           &
     & fv_cos_theta_latitude,                                           &
     & cf_area, cf_bulk, cf_liquid, cf_frozen,                          &
     & pstar, ntml, cumulus,                                            &
     & STASHwork,                                                       &
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
     &                  l_mixing_ratio, ICODE, CMESSAGE)

Use ac_diagnostics_mod

      IMPLICIT NONE

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



      INTEGER       INT18        ! Dummy variable for STASH_MAXLEN(18)
      INTEGER       P_FIELDDA,Q_LEVELSDA
      INTEGER model_levelsda             ! copy of model_levels
                                         !    for dynamic allocation
      INTEGER       ICODE        ! Return code : 0 Normal Exit
!                                !             : > 0 Error

      REAL theta_star_inc_halo(1-offx:row_length+offx,                  &
                                                       ! latest
     &               1-offy:rows+offy, model_levelsda) ! theta
      REAL q_star(P_FIELDDA,Q_LEVELSDA)            ! specific humidity,
      REAL qcl_star(P_FIELDDA,Q_LEVELSDA)          ! cloud liquid,
      REAL qcf_star(P_FIELDDA,Q_LEVELSDA)          ! cloud ice content
      CHARACTER*(*) CMESSAGE     ! Error message if return code >0
      INTEGER :: obs_flag_len,obs_len
      INTEGER :: OBS_FLAG(obs_flag_len)
      REAL    :: OBS(obs_len)
      Logical, intent(in)::                                             &
     & l_mixing_ratio            ! Use mixing ratio (if code available)

      Real, Intent (InOut) ::                                           &
     &  p(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &      model_levels+1)                                             &
     &, p_theta_levels(1-offx:row_length+offx,                          &
     &                   1-offy:rows+offy, model_levels)                &
     &, exner_theta_levels(1-offx:row_length+offx,                      &
     &                     1-offy:rows+offy, model_levels)              &
     &, r_theta_levels(1-offx:row_length+offx,                          &
     &                 1-offy:rows+offy, model_levels)                  &
     &, fv_cos_theta_latitude(1-offx:row_length+offx,                   &
     &                        1-offy:rows+offy)                         &
     &, cf_area(row_length, rows, wet_levels)                           &
     &, cf_bulk(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &           wet_levels)                                            &
     &, cf_liquid(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &             wet_levels)                                          &
     &, cf_frozen(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &             wet_levels)                                          &
     &, pstar(row_length,rows)                                          &
     &, ntml(row_length,rows)                                           &
     &, cumulus(row_length,rows) 

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
! CSIZEOBS start
! The variables involved in dimensioning observational data arrays
! in the data assimilation section are stored in this comdeck.

      ! For ATMOSPHERE assimilations the values are computed in the
      ! initialisation routine INITAC and then passed to the main
      ! routine AC.

      INTEGER :: A_MAX_NO_OBS    !  No of observations in AC Obs files.
      INTEGER :: A_MAX_OBS_SIZE  !  No of obs values in AC Obs files.

      COMMON /CSIZEOBS/ A_MAX_NO_OBS, A_MAX_OBS_SIZE
! CSIZEOBS end
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
! ACPARM start
! parameters used for dimensioning small permanent arrays
! values are maximum likely to avoid recompile on resolution change
! the actual dimensions are passed as arguments
! 5.0 17/06/99 Replace amaxsize by amaxsacp for expediency - to retain
!              old, invalid sizes temporarily. R. Rawlins
! ----------------------- include file: AMAXSACM -----------------------
! Description: Quick fix replacement of AMAXSIZE reference in ACPARM.
!              Contains sizes superseded at 5.0, but still embedded
!              throughout AC assimilation code. This allows AC
!              routines to compile without having to re-analyse and
!              correct every routine at this stage. This file should
!              only be of transient use as an interim solution.
!
! Current Code Owner: R. Rawlins
!         5.2   30/11/00 remove ROW_LENGTH_MAX and HORIZ_DIM_MAX
!                        already in amaxsize (now called in ac_ctl)
!                                          B Macpherson
!         6.1   31/08/04 Allow up to 100 levels.  R.Barnes
!         6.2   25/11/05 Set p_rows_max for N320.  R.Barnes

      INTEGER,PARAMETER:: P_ROWS_MAX = 481 ! Max number of rows
      INTEGER,PARAMETER:: P_LEVELS_MAX = 100 ! Max no. of total levels
      INTEGER,PARAMETER:: Q_LEVELS_MAX = 100 ! Max no. of wet levels

! AMAXSACP end
      INTEGER,PARAMETER:: NOBTYPMX = 182

! ACPARM end
!-----------------------------------------------------------------------
!*LCOMDECK COMACP
! Description:
!   Declares variables for the common block COMACP. This controls the
!  execution of the assimilation code.
!
! History:    see DOCACP
! 6.2      23/02/06     New logical REMOVE_NEG_LH. Mark Dixon
!-----------------------------------------------------------------------

! Declare parameters:
      INTEGER      MODEACP
        PARAMETER (MODEACP = 36)
      INTEGER      NANALTYP
        PARAMETER (NANALTYP = 30)
      INTEGER      NRADARS
        PARAMETER (NRADARS  = 15)

! Declare variables:
      INTEGER NACT,                        NPROG
      INTEGER AC_OBS_TYPES(NOBTYPMX),      LACT(NOBTYPMX)
      INTEGER GROUP_INDEX(NOBTYPMX),       TYPE_INDEX(NOBTYPMX)
      INTEGER GROUP_FIRST(NOBTYPMX),       GROUP_LAST(NOBTYPMX)
      INTEGER OBS_UNITNO,                  OBS_FORMAT
      INTEGER NO_OBS_FILES,                DIAG_RDOBS
      INTEGER IUNITNO,                     MGEOWT
      INTEGER N_GROUPS,                    GROUP_NO(NOBTYPMX)
      INTEGER MHCORFN,                     MACDIAG(MODEACP)
      INTEGER MWTFN,                       MDATADFN
      INTEGER NPASS_RF,                    NSLABS_SCFACT(P_LEVELS_MAX)
      INTEGER NO_SCFACT(NOBTYPMX),         IOMITOBS(NANALTYP)
      INTEGER MASTER_AC_TYPES(NOBTYPMX),   DEF_AC_ORDER(NOBTYPMX)
      INTEGER DEF_NO_ITERATIONS(NOBTYPMX), DEF_INTERVAL_ITER(NOBTYPMX)
      INTEGER DEF_NO_ANAL_LEVS(NOBTYPMX),  DEF_NO_WT_LEVS(NOBTYPMX)
      INTEGER DEF_MODE_HANAL(NOBTYPMX),    LENACT(NOBTYPMX)
      INTEGER DEF_OBTHIN(NOBTYPMX),        MVINT205
      INTEGER MRAMPFN,                     MGLOSSFN
      INTEGER      LHN_RANGE
      INTEGER      NPASS_RF_LHN


      REAL    OBTIME_NOM
      REAL    VERT_FILT
      REAL    GEOWT_H(P_ROWS_MAX -1)
      REAL    TROPLAT
      REAL    GEOWT_V(P_LEVELS_MAX)
      REAL    VERT_COR_SCALE(P_LEVELS_MAX, 4)
      REAL    VERT_CUTOFF_SL
      REAL    VERT_CUTOFF_BW
      REAL    VERT_CUTOFF_BH
      REAL    NON_DIV_COR
      REAL NON_DIV_COR_10M
      REAL    SPEED_LIMIT305
      REAL    TROPINT
      REAL    TIMEF_START
      REAL    TIMEF_OBTIME
      REAL    TIMEF_END
      REAL    CSCFACT_H(P_ROWS_MAX)
      REAL    CSCFACT_V(P_LEVELS_MAX)
      REAL    DEF_TIMEB(NOBTYPMX)
      REAL    DEF_TIMEA(NOBTYPMX)
      REAL    DEF_TGETOBB(NOBTYPMX)
      REAL    DEF_TGETOBA(NOBTYPMX)
      REAL    DEF_CSCALE_START(NOBTYPMX)
      REAL    DEF_CSCALE_OBTIME(NOBTYPMX)
      REAL    DEF_CSCALE_END(NOBTYPMX)
      REAL    DEF_RADINF(NOBTYPMX)
      REAL    WB_LAT_CC(P_ROWS_MAX)
      REAL    WB_VERT_V(P_LEVELS_MAX)
      REAL    WB_LAND_FACTOR
      REAL         RADAR_LAT(NRADARS)
      REAL         RADAR_LON(NRADARS)
      REAL         RADAR_RANGE_MAX
      REAL         EPSILON_LHN
      REAL         RELAX_CF_LHN
      REAL         F1_506 , F2_506 , F3_506
      REAL         ALPHA_LHN
      REAL         LHN_LIMIT
      REAL         FI_SCALE_LHN

      REAL    DEF_NUDGE_NH(NOBTYPMX)
      REAL    DEF_NUDGE_TR(NOBTYPMX)
      REAL    DEF_NUDGE_SH(NOBTYPMX)

      REAL    DEF_FI_VAR_FACTOR(NOBTYPMX)
      REAL    FI_SCALE
      REAL    FI_SCALE_FACTOR(P_LEVELS_MAX)
      REAL    DF_SCALE
      REAL    DF_SCALE_LEV(P_LEVELS_MAX)
      REAL    DF_COEFF(P_LEVELS_MAX)
      REAL    THRESH_DL
      REAL    THRESH_LM
      REAL    THRESH_MH
      REAL    THRESH_RMSF
      REAL    RADAR_RANGE
      REAL    NORTHLAT, SOUTHLAT, WESTLON, EASTLON
      REAL    VERT_COR_AERO

      LOGICAL LGEO
      LOGICAL LHYDR
      LOGICAL LHYDROL
      LOGICAL LSYN
      LOGICAL LTIMER_AC
      LOGICAL LAC_UARS
      LOGICAL LAC_MES
      LOGICAL LWBAL_SF,     LWBAL_UA
      LOGICAL WB_THETA_UA, WB_LAND_SCALE, WB_THETA_SF
      LOGICAL LRADAR (NRADARS)
      LOGICAL L_LATLON_PRVER
      LOGICAL L_MOPS_EQUALS_RH
      LOGICAL LCHECK_GRID
      LOGICAL      L_506_OBERR
      LOGICAL      L_LHN , L_LHN_SCALE
      LOGICAL      L_LHN_SEARCH , LHN_DIAG
      LOGICAL      L_VERIF_RANGE
      LOGICAL      L_LHN_LIMIT
      LOGICAL      L_LHN_FACT
      LOGICAL      L_LHN_FILT
      LOGICAL L_OBS_CHECK
      LOGICAL  REMOVE_NEG_LH
      LOGICAL USE_CONV_IN_MOPS

      COMMON /COMACP/ NACT,N_GROUPS,NPROG,                              &
     &  AC_OBS_TYPES,     LACT,              GROUP_NO,                  &
     &  LENACT,           LWBAL_SF,          LWBAL_UA,                  &
     &  LTIMER_AC,        LGEO,              LHYDR,                     &
     &  MGEOWT,           LSYN,              LAC_UARS,                  &
     &  OBS_UNITNO,       OBS_FORMAT,        NO_OBS_FILES,              &
     &  L_OBS_CHECK,                                                    &
     &  DIAG_RDOBS,       IUNITNO,           MVINT205,                  &
     &  MHCORFN,          MACDIAG,                                      &
     &  DEF_AC_ORDER,     DEF_NO_ITERATIONS, DEF_INTERVAL_ITER,         &
     &  MWTFN,            MDATADFN,          NSLABS_SCFACT,             &
     &  NO_SCFACT,        NPASS_RF,          MRAMPFN,                   &
     &  IOMITOBS,         TROPINT,           SPEED_LIMIT305,            &
     &  GEOWT_H,          GEOWT_V,           MGLOSSFN,                  &
     &  NON_DIV_COR,      TROPLAT,           VERT_FILT,                 &
     &  NON_DIV_COR_10M,                                                &
     &  VERT_COR_SCALE,                                                 &
     &  VERT_CUTOFF_SL,   VERT_CUTOFF_BW,    VERT_CUTOFF_BH,            &
     &  TIMEF_START,      TIMEF_OBTIME,      TIMEF_END,                 &
     &  CSCFACT_H,        CSCFACT_V,                                    &
     &  MASTER_AC_TYPES,                                                &
     &  DEF_NO_ANAL_LEVS, DEF_NO_WT_LEVS,    DEF_MODE_HANAL,            &
     &  DEF_TIMEB,        DEF_TIMEA,         DEF_TGETOBB,               &
     &  DEF_TGETOBA,      OBTIME_NOM,        DEF_OBTHIN,                &
     &  DEF_RADINF,       DEF_CSCALE_START,  DEF_CSCALE_OBTIME,         &
     &  DEF_CSCALE_END,                                                 &
     &  DEF_NUDGE_NH,     DEF_NUDGE_TR,      DEF_NUDGE_SH,              &
     &  GROUP_INDEX,      GROUP_FIRST,       GROUP_LAST,                &
     &  TYPE_INDEX,                                                     &
     &  FI_SCALE,         FI_SCALE_FACTOR,   DEF_FI_VAR_FACTOR,         &
     &  DF_SCALE,         DF_SCALE_LEV,      DF_COEFF,                  &
     &  LAC_MES,                                                        &
     &  THRESH_DL,        THRESH_LM,         THRESH_MH,                 &
     &  THRESH_RMSF,                                                    &
     &  RADAR_RANGE,      LRADAR,            LHYDROL,                   &
     &  L_LATLON_PRVER,   NORTHLAT,          SOUTHLAT,                  &
     &  WESTLON,          EASTLON,           L_MOPS_EQUALS_RH,          &
     &  LHN_RANGE ,  L_LHN , L_LHN_SCALE ,                              &
     &  L_LHN_SEARCH , LHN_DIAG , REMOVE_NEG_LH,                        &
     &  RADAR_LAT , RADAR_LON , RADAR_RANGE_MAX ,                       &
     &  EPSILON_LHN , ALPHA_LHN , RELAX_CF_LHN , LHN_LIMIT ,            &
     &  F1_506 , F2_506 , F3_506 ,                                      &
     &  L_506_OBERR , L_VERIF_RANGE , L_LHN_LIMIT , L_LHN_FACT ,        &
     &  L_LHN_FILT , FI_SCALE_LHN , NPASS_RF_LHN ,                      &
     &  VERT_COR_AERO,    LCHECK_GRID,                                  &
     &  WB_LAT_CC,        WB_VERT_V,         WB_LAND_FACTOR,            &
     & WB_THETA_UA,      WB_LAND_SCALE,   WB_THETA_SF, USE_CONV_IN_MOPS

! COMACP end

      REAL lcrcp
      PARAMETER (lcrcp=lc/cp)

! PC2 options are only important if L_pc2 = .true.
! Seek advice from the PC2 team before altering these parameters
! from .true., you will need to have put in place large amounts
! of extra code first.
      LOGICAL,PARAMETER:: L_pc2_cond=.true.  ! Do condensation and
! liquid cloud fraction changes.
      LOGICAL,PARAMETER:: L_pc2_cfl =.true.  ! Do liquid cloud
! fraction changes. This requires that the condensation as a
! result of assimilation is calculated directly.
! Note: One must not try to run with condensation on but liquid cloud
! fraction off with this code.
      LOGICAL,PARAMETER:: L_pc2_cff =.true.  ! Do ice cloud fraction
! changes. This requires that the ice increment from assimilation is
! calculated directly.

!L Dynamically allocated area for stash processing
      REAL STASHWORK(*)


      INTEGER                                                           &
     &       STASHMACRO_TAG,                                            &
                                       ! STASHmacro tag number
     &       MDI,                                                       &
                                       ! Missing data indicator
     &       K, ERROR                  ! do loop variable/ error

      INTEGER J                        ! DO Loop Variable.
      INTEGER I, IIND                  ! temporary scalars
      INTEGER im_index                 ! Internal model index
      INTEGER Q_LEVELS, P_FIELD          ! copies of Q_LEVELSDA,
                                         ! P_FIELDDA to avoid edits
      INTEGER rhc_row_length, rhc_rows   ! rhcrit dimensions
      INTEGER int_p                      ! field pointers
      INTEGER j_start, j_end             ! start and end indices
      INTEGER i_start, i_end             ! for processors
      INTEGER ji                         ! 2D array index for halo
                                         ! i/j variables
      INTEGER levels_per_level, large_levels ! needed by ls_arcld
      REAL RHCPT (1,1,Q_LEVELSDA)        ! rhcrit array
      REAL p_layer_centres(row_length, rows, 0:model_levelsda)
!          pressure at layer centres. Same as p_theta_levels
!          except bottom level = pstar, and at top = 0
      REAL p_layer_boundaries(row_length, rows, 0:model_levelsda)
!                pressure at layer boundaries. Same as p except at
!                bottom level = pstar, and at top = 0.
      REAL exner_layer_centres(row_length, rows, model_levelsda)
!          exner at layer centres. Same as exner_theta_levels
      REAL WORK(P_FIELDDA,Q_LEVELSDA)    ! array for large-scale
                                         ! latent heating
                                         ! or ls_cld dummy output
      REAL WORK2(P_FIELDDA,Q_LEVELSDA)   ! convective heating
                                         ! or ls_cld dummy output
      REAL theta_star(P_FIELDDA,model_levelsda) ! non-halo theta
      REAL bulk_cloud_nohalo(row_length,rows,q_levelsda)

      REAL DUMMY_FIELD(P_FIELDDA)

      ! These variables are for PC2 calculations
      REAL, DIMENSION(:,:,:), ALLOCATABLE::                             &
     &           q_work                                                 &
                                     ! Vapour work array (kg kg-1)
     &,          qcl_work                                               &
                                     ! Liquid work array (kg kg-1)
     &,          qcf_work                                               &
                                     ! Ice work array    (kg kg-1)
     &,          t_work                                                 &
                                     ! Temperature/theta work array
!                                    ! please see inline comments (K)
     &,          bulk_cloud_fraction                                    &
                                       ! Bulk cloud fraction (no units)
     &,          cloud_fraction_liquid                                  &
                                       ! Liquid cloud fraction  "
     &,          cloud_fraction_frozen                                  &
                                       ! Ice cloud fraction     "
     &,          delta_q                                                &
                                 ! Change in q to force condensation
     &,          delta_qcl                                              &
                                 ! Change in qcl to force condensation
     &,          delta_qcf                                              &
                                 ! Change in qcf
     &,          delta_t                                                &
                                 ! Change in t to force condensation
     &,          pc2_work        ! A work array for PC2.

! Declare allocatable arrays for passing cloud fractions
! to LS_ACF_Brooks
      Real, DIMENSION (:,:,:), ALLOCATABLE::                            &
     & cf_bulk_nohalo, cf_liquid_nohalo, cf_frozen_nohalo
!
! Start Routine
!
      Q_LEVELS = Q_LEVELSDA
      P_FIELD  = P_FIELDDA
      DO J=1,P_FIELDDA
         DUMMY_FIELD(J) =0.0
      END DO
!L
!L 1.0 Get address for each field from its STASH section/item code
!L     and STASHmacro tag  (searching only on STASHmacro tag)
      MDI            = IMDI
      STASHMACRO_TAG = 30

! Initialise STASHWORK for section 18.
      DO J = 1, INT18
        STASHWORK(J) = RMDI

      END DO

! Create local versions of exner and p on layer centres/boundaries.
! These local arrays have NO halo.
      Do j=1,rows
        Do i=1,row_length
          p_layer_centres(i,j,0) = pstar(i,j)
          p_layer_boundaries(i,j,0) = pstar(i,j)
        End Do
      End Do

      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i=1,  row_length
            p_layer_boundaries(i,j,k) = p(i,j,k+1) 
            p_layer_centres(i,j,k) = p_theta_levels(i,j,k) 
            exner_layer_centres(i,j,k)=exner_theta_levels(i,j,k) 
          End Do
        End Do
      End Do
      k=model_levels
      Do j = 1, rows
        Do i=1,  row_length
          p_layer_boundaries(i,j,k) = 0.0
          p_layer_centres(i,j,k) = p_theta_levels(i,j,k) 
          exner_layer_centres(i,j,k)=exner_theta_levels(i,j,k) 
        End Do
      End Do

!L 1.5 large scale rainfall rate LSRR stored in ac_diagnostics module

      IF (LSRR(1)  ==  rmdi) THEN
        ICODE    = 4203
        CMESSAGE = "AC_CTL: large scale rainfall rate not available"
        write(6,*)'AC_CTL 4203 ',ICODE,CMESSAGE
      END IF

      IF (ICODE  >   0) GOTO 999

!L 1.6 large scale snowfall rate LSSR stored in ac_diagnostics module

      IF (LSSR(1)  ==  rmdi) THEN
        ICODE    = 4204
        CMESSAGE = "AC_CTL: large scale snowfall rate not available"
        write(6,*)'AC_CTL 4204 ',ICODE,CMESSAGE
      END IF

      IF (ICODE  >   0) GOTO 999
      
      IF (USE_CONV_IN_MOPS) THEN 

!L 1.7 convective rainfall rate CVRR stored in ac_diagnostics module

       IF (CVRR(1)  ==  rmdi) THEN
         ICODE    = 5205
         CMESSAGE = "AC_CTL: convective rainfall rate not available"

       END IF

       IF (ICODE  >   0) GOTO 999

!L 1.8 convective snowfall rate CVSR stored in ac_diagnostics module

       IF (CVSR(1)  ==  rmdi) THEN
         ICODE    = 5206
         CMESSAGE = "AC_CTL: convective snowfall rate not available"

       END IF

       IF (ICODE  >   0) GOTO 999

!L 1.10 convective cloud cover on each model level CONVCC stored in module

       IF (CONVCC(1,1)  ==  rmdi) THEN
         ICODE    = 5212
         CMESSAGE = "AC_CTL: convective cloud amount not available"

       END IF

       IF (ICODE  >   0) GOTO 999

!L 1.11 bulk cloud fraction (liquid+ice) after large scale cloud CF_LSC 

       IF (CF_LSC(1,1)  ==  rmdi) THEN
         ICODE    = 9201
         CMESSAGE = "AC_CTL: cf after large scale cloud not available"
       ENDIF

       IF (ICODE  >   0) GOTO 999

      ENDIF  !(IF USE_CONV_IN_MOPS)
     
      IF( L_LHN ) THEN
!  seek convective heating rate and
!  diagnostics for calculating large-scale latent heating rate

      IF (USE_CONV_IN_MOPS) THEN 
!L 1.13 temperature increment across convection TINC_CVN stored in module

       IF (TINC_CVN(1,1)  ==  rmdi) THEN
         ICODE    = 5181
         CMESSAGE = "AC_CTL: temp incrs across conv'n not available"
       ENDIF

       IF (ICODE  >   0) GOTO 999

      ENDIF  

!L 1.14 temperature increment across large scale precipitation TINC_PPN

      IF (TINC_PPN(1,1)  ==  rmdi) THEN
        ICODE    = 4181
        CMESSAGE = "AC_CTL: temp incrs across ls_ppn not available"
      ENDIF

      IF (ICODE  >   0) GOTO 999

!L 1.15 cloud liquid water after large scale cloud QCL_LSC stored in module

      IF (QCL_LSC(1,1)  ==  rmdi) THEN
        ICODE    = 9206
        CMESSAGE = "AC_CTL: qcl after large scale cloud not available"
      ENDIF

      IF (ICODE  >   0) GOTO 999

!L 1.16 cloud liquid water after advection QCL_ADV stored in ac_diagnostics module

      IF (QCL_ADV(1,1)  ==  rmdi) THEN
        ICODE    = 12254
        CMESSAGE = "AC_CTL: qcl after advection not available"
      ENDIF

      IF (ICODE  >   0) GOTO 999


! 1.25  Calculate 'large-scale' latent heating contributions
!       ----------------------------------------------------
       DO K=1,Q_LEVELS
         DO J=1,P_FIELD
          WORK(J,K) = TINC_PPN(J,K) +                                   &
     &              lcrcp*( QCL_LSC(J,K)- QCL_ADV(J,K)  )
         END DO
       END DO
!  large scale latent heating currently dT/dt in K/timestep
!  as is convective heating - convert both to dtheta/dt in K/s
       IF (USE_CONV_IN_MOPS) THEN
        Do k = 1, q_levels
          Do j = 1, rows
            Do i=1,  row_length
              int_p = (j-1) * row_length + i
              WORK(int_p,k)  =  WORK(int_p,k) /                          &
     &           (exner_layer_centres(i,j,k) * SECS_PER_STEPim(atmos_im))
              WORK2(int_p,k) =  TINC_CVN(int_p,k) /                      &
     &           (exner_layer_centres(i,j,k) * SECS_PER_STEPim(atmos_im))
            End Do
          End Do
        End Do

       ELSE     ! IF USE_CONV_IN_MOP is false

        Do k = 1, q_levels
          Do j = 1, rows
            Do i=1,  row_length
              int_p = (j-1) * row_length + i
              WORK(int_p,k)  =  WORK(int_p,k) /                          &
     &           (exner_layer_centres(i,j,k) * SECS_PER_STEPim(atmos_im))
              WORK2(int_p,k) = 0.0
            End Do
          End Do
        End Do
       ENDIF

      ELSE     !  if LHN not selected
!                 initialise dummy heating rate array to pass to AC
       DO K=1,Q_LEVELS
         DO J=1,P_FIELD
           WORK(J,K) =  0.0
           WORK2(J,K) =  0.0
         END DO
       END DO

      END IF   !  L_LHN

!L----------------------------------------------------------------------
!L 2. --- Section 18 Data Assimilation ------

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('AC      ', 3)

      END IF

      im_index=internal_model_index(A_IM)


!  copy non halo values from theta_star_inc_halo to
!  2-d array theta_star
      Do K = 1, model_levels
        Do J = 1, rows
          DO I=1,  row_length
            int_p = (J-1) * row_length + I
            theta_star(int_p,K) = theta_star_inc_halo(I,J,K)
          End Do
        End Do
      End Do

        If (L_pc2 .and.                                                 &
     &     (L_pc2_cond .or. L_pc2_cfl .or. L_pc2_cff)   ) then
! Input values needed for the PC2 scheme
        Allocate ( q_work  (row_length,rows,q_levelsda) )
        Allocate ( qcl_work(row_length,rows,q_levelsda) )
        Allocate ( qcf_work(row_length,rows,q_levelsda) )
        Allocate ( t_work  (row_length,rows,q_levelsda) )
        Allocate ( pc2_work(row_length,rows,q_levelsda) )
        Do k = 1, q_levelsda
          Do j = 1, rows
            Do i= 1, row_length
              int_p = (j-1) * row_length + i
              ! Copy input values into the work arrays
              q_work  (i,j,k) = q_star  (int_p,k)
              qcl_work(i,j,k) = qcl_star(int_p,k)
              qcf_work(i,j,k) = qcf_star(int_p,k)
              t_work  (i,j,k) = theta_star(int_p,k)
              ! t_work is now theta before AC assimilation
              pc2_work(i,j,k) = exner_layer_centres(i,j,k)
              ! pc2_work is now exner before AC assimilation
            End Do
          End Do
        End Do
      End If  ! L_pc2 .and. (L_pc2_cond.or.L_pc2_cfl.or.L_pc2_cff)

      IF (USE_CONV_IN_MOPS) THEN

! DEPENDS ON: ac
       CALL AC (                                                         &
     &   model_levels, wet_levels, row_length, rows, BL_LEVELS,          &
     &   A_MAX_OBS_SIZE, A_MAX_NO_OBS, theta_field_size,                 &
     &   STEPim(atmos_im), SECS_PER_STEPim(atmos_im),                    &
     &   exner_layer_centres, pstar,                                     &
     &   p_layer_centres(1,1,1),                                         &
     &   theta_star, q_star, qcl_star, qcf_star,                         &
     &   CONVCC, LSRR, LSSR, CVRR, CVSR,                                 &
     &   CF_LSC,WORK2,WORK, RHCRIT, L_eacf,                              &
     &   OBS_FLAG,OBS,                                                   &
     &   STINDEX(1,1,18,im_index),                                       &
     &   STLIST, LEN_STLIST, SI(1,18,im_index), SF(1,18),                &
     &   STASHWORK, STASH_LEVELS,                                        &
     &   NUM_STASH_LEVELS, STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,        &
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
     &   ICODE, CMESSAGE)
      ELSE
       CALL AC (                                                         &
     &   model_levels, wet_levels, row_length, rows, BL_LEVELS,          &
     &   A_MAX_OBS_SIZE, A_MAX_NO_OBS, theta_field_size,                 &
     &   STEPim(atmos_im), SECS_PER_STEPim(atmos_im),                    &
     &   exner_layer_centres, pstar,                                     &
     &   p_layer_centres(1,1,1),                                         &
     &   theta_star, q_star, qcl_star, qcf_star,                         &
     &   WORK2, LSRR, LSSR, DUMMY_FIELD, DUMMY_FIELD,                                 &
     &   CF_LSC,WORK2,WORK, RHCRIT, L_eacf,                              &
     &   OBS_FLAG,OBS,                                                   &
     &   STINDEX(1,1,18,im_index),                                       &
     &   STLIST, LEN_STLIST, SI(1,18,im_index), SF(1,18),                &
     &   STASHWORK, STASH_LEVELS,                                        &
     &   NUM_STASH_LEVELS, STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,        &
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
     &   ICODE, CMESSAGE)
      ENDIF

      ! For PC2 we need to calculate the forcing values of qT and TL
        If (L_pc2) then

        If (L_pc2_cond .or. L_pc2_cfl .or. L_pc2_cff) then

        Allocate ( delta_q  (row_length,rows,q_levelsda) )
        Allocate ( delta_qcl(row_length,rows,q_levelsda) )
        Allocate ( delta_qcf(row_length,rows,q_levelsda) )
        Allocate ( delta_t  (row_length,rows,q_levelsda) )
        Allocate ( bulk_cloud_fraction(row_length,rows,q_levelsda)  )
        Allocate ( cloud_fraction_liquid(row_length,rows,q_levelsda))
        Allocate ( cloud_fraction_frozen(row_length,rows,q_levelsda))
        Do k = 1, q_levelsda
          Do j = 1, rows
            Do i= 1, row_length

              int_p = (j-1) * row_length + i
              ! Calculate the increments across the AC scheme
              delta_q  (i,j,k) =  q_star  (int_p,k) - q_work  (i,j,k)
              delta_qcl(i,j,k) =  qcl_star(int_p,k) - qcl_work(i,j,k)
              delta_qcf(i,j,k) =  qcf_star(int_p,k) - qcf_work(i,j,k)

              ! t_work is currently theta before AC assimilation (K)
              delta_t  (i,j,k) =  exner_layer_centres(i,j,k)            &
     &                         * (theta_star(int_p,k) - t_work(i,j,k) )
              ! deltat is now change in temperature associated with
              ! AC assimilation minus the adiabatic contribution
              ! associated with the change in pressure. This last
              ! part is dealt with in pc2_pressure_ctl.
              ! Pc2_work currently contains exner before AC assim
              t_work  (i,j,k) =  pc2_work(i,j,k) * t_work(i,j,k)
              ! t_work now contains temperature before AC assimilation

              pc2_work(i,j,k) =  0.0
              ! pc2_work now contains an array of zeros

              ! Now copy cloud fraction information from the prognostics
              ! arrays
              bulk_cloud_fraction  (i,j,k) = cf_bulk(i,j,k) 
              cloud_fraction_liquid(i,j,k) = cf_liquid(i,j,k) 
              cloud_fraction_frozen(i,j,k) = cf_frozen(i,j,k) 

            End Do
          End Do
        End Do

        ! Now force condensation and cloud fraction updates.
        ! PC2 is not set up for prognostic area_cloud_fraction.

! DEPENDS ON: pc2_assim
        Call pc2_assim( q_levelsda, row_length, rows                    &
     &,                 SECS_PER_STEPim(atmos_im), l_pc2_cfl, l_pc2_cff &
     &,                 l_mixing_ratio                                  &
     &,                 t_work, bulk_cloud_fraction                     &
     &,                 cloud_fraction_liquid, cloud_fraction_frozen    &
     &,                 q_work, qcl_work, qcf_work                      &
     &,              p_layer_centres(1:row_length,1:rows,1:q_levelsda)  &
     &,                 delta_t, delta_q                                &
     &,                 delta_qcl, delta_qcf, pc2_work                  &
     &                  )

        ! q_work, qcl_work and t_work have now been updated by the
        ! forcing and condensation terms together.
        ! cloud_fraction_liquid (and _frozen and _bulk) has been
        ! updated by the condensation. qcf_work and p_work are not
        ! updated.

        ! Now copy q_work, qcl_work and t_work variables back
        ! into the inout variables if we haven't a confident
        ! direct estimate of condensation from another source.
        If (L_pc2_cond) then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i= 1, row_length
                int_p = (j-1) * row_length + i
                q_star  (int_p,k) = q_work  (i,j,k)
                qcl_star(int_p,k) = qcl_work(i,j,k)
                ! Remember the inout variable is theta, not temp.
                theta_star(int_p,k) = t_work(i,j,k)                     &
     &                              / exner_layer_centres(i,j,k)
              End Do
            End Do
          End Do
        End If  ! L_pc2_cond

        ! Update the D1 cloud fractions
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i= 1, row_length
              cf_bulk(i,j,k) = bulk_cloud_fraction(i,j,k) 
              cf_liquid(i,j,k) = cloud_fraction_liquid(i,j,k) 
              cf_frozen(i,j,k) = cloud_fraction_frozen(i,j,k) 
            End Do
          End Do
        End Do

          If (L_CLD_AREA) Then
            If (L_ACF_Brooks) Then
              Allocate ( cf_bulk_nohalo(row_length,rows,wet_levels) )
              Allocate ( cf_liquid_nohalo(row_length,rows,wet_levels) )
              Allocate ( cf_frozen_nohalo(row_length,rows,wet_levels) )

              i_start = 1
              i_end = row_length
              j_start = 1
              j_end = rows
              If (model_domain  ==  mt_lam) Then
                If(at_extremity(PSouth)) j_start = 2
                If(at_extremity(PNorth)) j_end = rows-1
                If(at_extremity(PWest)) i_start = 2
                If(at_extremity(PEast)) i_end = row_length-1
              End If
              If (model_domain  ==  mt_cyclic_lam) Then
                If (at_extremity(PSouth)) j_start = 2
                If (at_extremity(PNorth)) j_end = rows-1
              End If

! Place bulk, liquid and frozen cloud fractions in halo-free arrays
              Do k = 1, wet_levels
                Do j = j_start, j_end
                  Do i = i_start, i_end
                    ji = i+halo_i + (j+halo_j-1) * (row_length+2*halo_i)
                    cf_bulk_nohalo(i,j,k)  = cf_bulk(i,j,k)
                    cf_liquid_nohalo(i,j,k)= cf_liquid(i,j,k)
                    cf_frozen_nohalo(i,j,k)= cf_frozen(i,j,k)
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

        ! Now deallocate the arrays
        Deallocate(q_work  )
        Deallocate(qcl_work)
        Deallocate(qcf_work)
        Deallocate(t_work  )
        Deallocate(pc2_work)
        Deallocate(delta_q  )
        Deallocate(delta_qcl)
        Deallocate(delta_qcf)
        Deallocate(delta_t  )
        Deallocate(bulk_cloud_fraction  )
        Deallocate(cloud_fraction_liquid)
        Deallocate(cloud_fraction_frozen)

        End If !  L_pc2_cond.or.L_pc2_cfl.or.L_pc2_cff

      Else  ! L_pc2

!  2.1 call cloud scheme to 'rebalance' thermodynamic fields
!  calculate Tl and qt
         Do k=1,q_levels
           Do j = 1, rows
             Do i=1,  row_length
               int_p = (j-1) * row_length + i
               theta_star(int_p,k) = theta_star(int_p,k)*               &
     &             exner_layer_centres(i,j,k) - lcrcp*qcl_star(int_p,k)
               q_star(int_p,k) = q_star(int_p,k) + qcl_star(int_p,k)
             End Do
           End Do
         End Do
! set up some arguments for ls_cld

      rhc_row_length = 1
      rhc_rows       = 1

        Do k = 1, q_levelsda
          Do j = 1, rows
            Do i= 1, row_length
              bulk_cloud_nohalo(i,j,k) = cf_bulk(i,j,k) 
            End Do
          End Do
        End Do
      DO K = 1, Q_LEVELS
        RHCPT (1,1,K) = RHCRIT(K)
      END DO
! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
        levels_per_level = 3
        large_levels = ((q_levels - 2)*levels_per_level) + 2

! DEPENDS ON: ls_arcld
        CALL ls_arcld( p_layer_centres, RHCPT, p_layer_boundaries,      &
     &                 model_levels, q_levels, row_length, rows,        &
     &                 rhc_row_length, rhc_rows, bl_levels,             &
     &                 cloud_fraction_method,overlap_ice_liquid,        &
     &                 ice_fraction_method,ctt_weight,t_weight,         &
     &                 qsat_fixed,sub_cld,                              &
     &                 levels_per_level, large_levels,                  &
     &                 L_cld_area,L_ACF_Cusack,L_ACF_Brooks,L_eacf,     &
     &                 halo_i, halo_j, offx, offy,                      &
     &                 delta_lambda, delta_phi,                         &
     &                 r_theta_levels, FV_cos_theta_latitude,           &
     &                 ntml, cumulus, l_mixing_ratio,                   &
     &                 qcf_star, theta_star, q_star, qcl_star,          &
     &                 CF_AREA(:,:,1), bulk_cloud_nohalo,               &
     &                 WORK2, WORK,                                     &
     &                 ERROR,mype)

        Do k = 1, wet_levels
          Do j = 1, rows
            Do i= 1, row_length
              cf_bulk(i,j,k) = bulk_cloud_nohalo(i,j,k) 
            End Do
          End Do
        End Do

! "1.24.4"  convert t back to theta
         Do k=1,q_levels
           Do j = 1, rows
             Do i=1,  row_length
               int_p = (j-1) * row_length + i
               theta_star(int_p,k) = theta_star(int_p,k)/               &
     &                         exner_layer_centres(i,j,k)
             End Do
           End Do
         End Do

      End If  ! L_pc2

!  copy back theta_star into theta_star_inc_halo
      Do K = 1, model_levels
        Do J = 1, rows
          DO I=1,  row_length
            int_p = (J-1) * row_length + I
            theta_star_inc_halo(I,J,K) = theta_star(int_p,K)
          End Do
        End Do
      End Do


        IF (LTIMER) THEN
! DEPENDS ON: timer
          CALL TIMER('AC      ', 4)

        END IF

!! DEPENDS ON: stash
!        CALL STASH(a_sm, a_im, 18, STASHWORK,                           &
!#include "argdumga.h"
!#include "argsts.h"
!#include "argppx.h"
!     &                                      ICODE, CMESSAGE)
!
!        IF (ICODE >  0) GOTO 999

 999  CONTINUE
      RETURN
      END SUBROUTINE AC_CTL
