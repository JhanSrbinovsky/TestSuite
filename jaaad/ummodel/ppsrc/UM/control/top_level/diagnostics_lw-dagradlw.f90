
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_lw
!***********************************************************************
!*  This section has the following development path                   **
!*  vn2.6  ...  Basic code                                            **
!*  vn2.6a ...  code fixed/developed                                  **
!***********************************************************************

      Subroutine diagnostics_lw(                                        &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels, ozone_levels             &
     &,                      cloud_levels                               &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      timestep                                   &
     &,                      T_n, T_inc                                 &
     &,                      q_n, qcl_n, cf_n, cfl_n                    &
     &,                      T_latest, q_latest, qcl_latest             &
     &,                      cf_latest, cfl_latest                      &
     &,                      surflw, OLR, clear_olr                     &
     &,                      total_cloud_cover_lw                       &
     &,                      surface_down_flux_lw, surf_down_clr_lw     &
     &,                      T_incr_diagnostic                          &
     &,                      clear_hr_lw                                &
     &,                      net_flux_trop_lw, down_flux_trop_lw        &
     &,                      total_cloud_on_levels                      &
     &,                      cloud_absorptivity                         &
     &,                      cloud_weight_absorptivity                  &
     &,                      ls_cloud_absorptivity                      &
     &,                      ls_cloud_weight_absorptivity               &
     &,                      cnv_cloud_absorptivity                     &
     &,                      cnv_cloud_weight_absorptivity              &
     &,                      isccp_weights                              &
     &,                      isccp_cf                                   &
     &,                      isccp_cf_tau_0_to_p3                       &
     &,                      isccp_cf_tau_p3_to_1p3                     &
     &,                      isccp_cf_tau_1p3_to_3p6                    &
     &,                      isccp_cf_tau_3p6_to_9p4                    &
     &,                      isccp_cf_tau_9p4_to_23                     &
     &,                      isccp_cf_tau_23_to_60                      &
     &,                      isccp_cf_tau_ge_60                         &
     &,                      meanalbedocld                              &
     &,                      meantaucld                                 &
     &,                      meanptop                                   &
     &,                      totalcldarea                               &
     &,                      ls_qcl_rad                                 &
     &,                      ls_qcf_rad                                 &
     &,                      cc_qcl_rad                                 &
     &,                      cc_qcf_rad                                 &
     &,                      ls_cl_rad                                  &
     &,                      ls_cf_rad                                  &
     &,                      cc_cl_rad                                  &
     &,                      cc_cf_rad                                  &
     &,                      ozone                                      &
     &,                      O3_trop_level                              &
     &,                      O3_trop_height                             &
     &,                      T_trop_level                               &
     &,                      T_trop_height                              &
     &,                      LW_incs, LWsea                             &
     &,                      n_aod_wavel                                &
     &,                      aod_sulphate                               &
     &,                      aod_dust                                   &
     &,                      aod_seasalt                                &
     &,                      aod_soot                                   &
     &,                      aod_biomass                                &
     &,                      aod_biogenic                               &
     &,                      aod_ocff                                   &
     &,                      aod_delta                                  &
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork                                                        &
     &     )


! Purpose:
!          Calculates diagnostics and outputs them.
!
! Method:
!
! History:
! Date     Version     Comment
! ----     -------     -------
! 22/04/98   2.4   Distinguish between area and bulk cloud fractions in
!                  radiation calculations: bulk for Re, area otherwise.
!                  Add Total Cloud in LW diagnostic. Andrew C. Bushell
! Version   Date     Comment
! ---      -----     -------
! 5.0  30/11/99 Original version in UM. J-C Thil.
! 5.1  09/12/99 Add error trapping. Rick Rawlins
! 5.1  28/02/00 Provide explicit model increments as STASH output
!               diagnostics. R Rawlins
! 5.2  14/11/00 Provide code for reactivated diagnostics.
!                                 (J. M. Edwards)
! 5.4  03/08/01 Diagnostics for PC2 cloud scheme. D.R. Wilson
! 5.5  09/01/02 Add absorptivity and isccp diagnostics.
!                                             A.B.Keen/K.D.Williams
!
! 5.5  24/02/03 Output O3_trop_level, O3_trop_height,
!               T_trop_level, T_trop_height.
!                                      (J.-C. Thelen)
! 6.1  05/05/04 Add extra ISCCP diagnostics.     A. Jones
! 6.2  05/11/05 Add aerosol optical depth diagnostics.    N. Bellouin
! 6.4  16/01/07 Add new grid-box mean cloud diagnostics.  Adrian Lock
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
      Real                                                              &
     &  timestep         ! atmosphere timestep


      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, cloud_levels                                                    &
                         ! number of cloudy levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, ozone_levels                                                    &
                         ! number of levels where ozone is held
     &, number_format                                                   &
                         ! switch controlling number format diagnostics
                         ! are written out in. See PP_WRITE for details.
     &, model_domain     ! indicator as to model type, ie global, lam

      Integer                                                           &
     &  global_row_length                                               &
                            !IN. NUMBER OF points on a global row
     &, global_rows                                                     &
                            !IN. NUMBER OF global rows
     &, me                                                              &
                            !IN. Processor number
     &, halo_i                                                          &
                            !IN. size of large halo in x direction
     &, halo_j                                                          &
                            !IN. size of large halo in y direction
     &, off_x                                                           &
                            !IN. size of small halo in x direction
     &, off_y                                                           &
                            !IN. size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)

! Primary Arrays used in all models
      Real                                                              &
     &  T_n(row_length, rows, model_levels)                             &
     &, T_inc(row_length, rows, model_levels)                           &
     &, q_n(row_length, rows, wet_model_levels)                         &
     &, qcl_n(row_length, rows, wet_model_levels)                       &
     &, cf_n(row_length, rows, wet_model_levels)                        &
     &, cfl_n(row_length, rows, wet_model_levels)                       &
     &, T_latest(row_length, rows, model_levels)                        &
     &, q_latest(row_length, rows, wet_model_levels)                    &
     &, qcl_latest(row_length, rows, wet_model_levels)                  &
     &, cf_latest(row_length, rows, wet_model_levels)                   &
     &, cfl_latest(row_length, rows, wet_model_levels)                  &
     &, LW_incs(row_length, rows, 0:model_levels)

      Real                                                              &
     &  OLR (row_length, rows)                                          &
     &, surflw (row_length, rows)                                       &
     &, clear_olr (row_length, rows)                                    &
     &, total_cloud_cover_lw (row_length, rows)                         &
     &, surface_down_flux_lw (row_length, rows)                         &
     &, surf_down_clr_lw (row_length, rows)                             &
     &, LWsea(row_length, rows)                                         &
     &, T_incr_diagnostic(row_length,rows,model_levels)                 &
     &, clear_hr_lw(row_length,rows,model_levels)                       &
     &, net_flux_trop_lw(row_length,rows)                               &
     &, down_flux_trop_lw(row_length,rows)                              &
     &, ozone(row_length,rows,ozone_levels)                             &
     &, total_cloud_on_levels(row_length, rows, cloud_levels)           &
     &, cloud_absorptivity(row_length, rows, cloud_levels)              &
     &, cloud_weight_absorptivity(row_length, rows, cloud_levels)       &
     &, ls_cloud_absorptivity(row_length, rows, cloud_levels)           &
     &, ls_cloud_weight_absorptivity(row_length, rows, cloud_levels)    &
     &, cnv_cloud_absorptivity(row_length, rows, cloud_levels)          &
     &, cnv_cloud_weight_absorptivity(row_length, rows, cloud_levels)   &
     &, isccp_weights(row_length, rows)                                 &
     &, isccp_cf(row_length, rows, 7)                                   &
     &, isccp_cf_tau_0_to_p3(row_length, rows, 7)                       &
     &, isccp_cf_tau_p3_to_1p3(row_length, rows, 7)                     &
     &, isccp_cf_tau_1p3_to_3p6(row_length, rows, 7)                    &
     &, isccp_cf_tau_3p6_to_9p4(row_length, rows, 7)                    &
     &, isccp_cf_tau_9p4_to_23(row_length, rows, 7)                     &
     &, isccp_cf_tau_23_to_60(row_length, rows, 7)                      &
     &, isccp_cf_tau_ge_60(row_length, rows, 7)                         &
     &, meanalbedocld(row_length,rows)                                  &
     &, meantaucld(row_length,rows)                                     &
     &, meanptop(row_length,rows)                                       &
     &, totalcldarea(row_length,rows)                                   &
     &, ls_qcl_rad(row_length,rows,model_levels)                        &
     &, ls_qcf_rad(row_length,rows,model_levels)                        &
     &, cc_qcl_rad(row_length,rows,model_levels)                        &
     &, cc_qcf_rad(row_length,rows,model_levels)                        &
     &, ls_cl_rad(row_length,rows,model_levels)                         &
     &, ls_cf_rad(row_length,rows,model_levels)                         &
     &, cc_cl_rad(row_length,rows,model_levels)                         &
     &, cc_cf_rad(row_length,rows,model_levels)                         &
     &, O3_trop_level(row_length,rows)                                  &
     &, O3_trop_height(row_length,rows)                                 &
     &, T_trop_level(row_length,rows)                                   &
     &, T_trop_height(row_length,rows)

! Aerosol optical depth diagnostics
      Integer                                                           &
     &  n_aod_wavel
      Real                                                              &
     &  aod_sulphate(row_length, rows, n_aod_wavel)                     &
     &, aod_dust(row_length, rows, n_aod_wavel)                         &
     &, aod_seasalt(row_length, rows, n_aod_wavel)                      &
     &, aod_soot(row_length, rows, n_aod_wavel)                         &
     &, aod_biomass(row_length, rows, n_aod_wavel)                      &
     &, aod_biogenic(row_length, rows, n_aod_wavel)                     &
     &, aod_ocff(row_length, rows, n_aod_wavel)                         &
     &, aod_delta(row_length, rows, n_aod_wavel)

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

! Diagnostic variables
      Real                                                              &
     & STASHwork(*)    ! STASH workspace

! Local array & variables
      Real                                                              &
     &  work_3d(row_length, rows, model_levels)                         &
     &  ,isccp_dummy_3d(row_length,rows,cloud_levels)

      Integer                                                           &
     &  i, j, k                                                         &
     &,    icode                                                        &
                                ! Return code  =0 Normal exit  >1 Error
     & ,item                                                            &
                        ! STASH item
     & ,sect            ! STASH section
      Parameter( sect = 2 ) ! for lw radiation

      Character*80  cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_lw')

      Integer                                                           &
     &  im_index        ! internal model index

! External routines
      External                                                          &
     &  copydiag, copydiag_3d                                           &
     &  ,Ereport

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------

! pseudo temperature after lw radiation (diagnostic only)

      item =   4  ! T + LW T_increment
      If (icode <= 0 .and. sf(item,sect)) Then

        Do k = 1,model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3d(i,j,k) = T_n(i,j,k) + T_inc(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 004)"//cmessage
        End if

      Endif  !  sf(item,sect)

! increment diagnostics= modified - previous

      item = 161  ! temperature increment
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        T_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 161)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 181  ! temperature increment including condensation
      If (icode <= 0 .and. sf(item,sect)) Then
!
      Do k = 1,model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = T_latest(i,j,k) - T_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 181)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 182  ! Vapour increment
      If (icode <= 0 .and. sf(item,sect)) Then
!
      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = q_latest(i,j,k) - q_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 182)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 183  ! liquid water content increment
      If (icode <= 0 .and. sf(item,sect)) Then
!
      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = qcl_latest(i,j,k) - qcl_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 183)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 192  ! total cloud fraction increment
      If (icode <= 0 .and. sf(item,sect)) Then
!
      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = cf_latest(i,j,k) - cf_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 192)"//cmessage
         End if

      Endif  !  sf(item,sect)
!
      item = 193  ! liquid cloud fraction increment
      If (icode <= 0 .and. sf(item,sect)) Then
!
      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = cfl_latest(i,j,k) - cfl_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 193)"//cmessage
         End if

      Endif  !  sf(item,sect)

      If (sf(201,2)) Then

!L   surflw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(201,2,im_index)),surflw,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,201,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if



      If (sf(205,2)) Then

!L   OLR

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(205,2,im_index)),OLR,              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,205,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 205)"
            goto 9999
         End if

      End if


      If (sf(206,2)) Then

!L  clear_olr

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(206,2,im_index)),clear_olr,        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,206,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 206)"
            goto 9999
         End if

      End if




      If (sf(203,2)) Then

!L    LWsea : 'net down lw rad flux: open sea'

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(203,2,im_index)),                  &
     &        LWsea,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,203,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 203 = LWsea)"
            goto 9999
         End if

      End if


      If (sf(204,2)) Then

!L  total_cloud_cover_lw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(204,2,im_index)),                  &
     &        total_cloud_cover_lw,                                     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,204,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 204)"
            goto 9999
         End if

      End if



      item = 207 ! surface downward LW flux
      If (icode <= 0 .and. sf(item,sect)) Then

!L    surface_down_flux_lw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(item,sect,im_index)),              &
     &        surface_down_flux_lw,                                     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 207)"
            goto 9999
         End if

      End if



      item = 208 ! surface downward clear skies LW flux
      If (icode <= 0 .and. sf(item,sect)) Then

!L    surf_down_clr_lw,

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(item,sect,im_index)),              &
     &        surf_down_clr_lw,                                         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
         End if

      End if




      If (sf(233,2)) Then

!L  clear_hr_lw

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(233,2,im_index)),               &
     &        clear_hr_lw,                                              &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,233,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,233,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 233)"//cmessage
            goto 9999
         End if

      End if



      If (sf(237,2)) Then

!L  net_flux_trop_lw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(237,2,im_index)),                  &
     &        net_flux_trop_lw,                                         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,237,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 237)"
            goto 9999
         End if

      End if



      If (sf(238,2)) Then

!L  down_flux_trop_lw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(238,2,im_index)),                  &
     &        down_flux_trop_lw,                                        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,238,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 238)"
            goto 9999
         End if

      End if


      If (sf(260,2)) Then

!L  Ozone

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(260,2,im_index)),               &
     &        ozone,                                                    &
     &        row_length,rows,ozone_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,260,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,260,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 260)"//cmessage
            goto 9999
         End if

      End if

      If (sf(280,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(280,2,im_index)),                  &
     &        O3_trop_level,                                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,280,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 280)"//cmessage
         End if
      End if

      If (sf(281,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(281,2,im_index)),                  &
     &        O3_trop_height,                                           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,281,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 281)"//cmessage
         End if
      End if

      If (sf(282,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(282,2,im_index)),                  &
     &        T_trop_level,                                             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,282,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 282)"//cmessage
         End if
      End if

      If (sf(283,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(283,2,im_index)),                  &
     &        T_trop_height,                                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,283,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 283)"//cmessage
         End if
      End if
! LW heating =
!            lw radiation temperature increment per timestep / timestep

      item = 232  ! LW heating
      If (icode <= 0 .and. sf(item,sect)) Then

        Do k = 1,model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3d(i,j,k) = LW_incs(i,j,k)/timestep
            Enddo ! i
          Enddo   ! j
        Enddo     ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 232)"//cmessage
        End if

      Endif  !  sf(item,sect)


      If (sf(261,2)) Then

!   total_cloud_on_levels

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(261,2,im_index)),               &
     &        total_cloud_on_levels,                                    &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,261,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,261,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 261)"//cmessage
            goto 9999
         End if

      End if

       If (sf(262,2)) Then

!   cloud_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(262,2,im_index)),               &
     &        cloud_absorptivity,                                       &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,262,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,262,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 262)"//cmessage
            goto 9999
         End if

      End if

      If (sf(263,2)) Then

!   cloud_weight_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(263,2,im_index)),               &
     &        cloud_weight_absorptivity,                                &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,263,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,263,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 263)"//cmessage
            goto 9999
         End if

      End if

      If (sf(264,2)) Then

!   ls_cloud_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(264,2,im_index)),               &
     &        ls_cloud_absorptivity,                                    &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,264,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,264,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 264)"//cmessage
            goto 9999
         End if

      End if

      If (sf(265,2)) Then

!   ls_cloud_weight_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(265,2,im_index)),               &
     &        ls_cloud_weight_absorptivity,                             &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,265,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,265,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 265)"//cmessage
            goto 9999
         End if

      End if

      If (sf(266,2)) Then

!   cnv_cloud_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(266,2,im_index)),               &
     &        cnv_cloud_absorptivity,                                   &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,266,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,266,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 266)"//cmessage
            goto 9999
         End if

      End if

      If (sf(267,2)) Then

!   cnv_cloud_weight_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(267,2,im_index)),               &
     &        cnv_cloud_weight_absorptivity,                            &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,267,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,267,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 267)"//cmessage
            goto 9999
         End if

      End if

      If (sf(269,2)) Then

!   isccp_weights

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(269,2,im_index)),                  &
     &        isccp_weights,                                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,269,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 269)"
            goto 9999
         End if

      End if

      If (sf(290,2)) Then

!   meanalbedocld

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(290,2,im_index)),                  &
     &        meanalbedocld,                                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,290,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 290)"
            goto 9999
         End if

      End if

      If (sf(291,2)) Then

!   meantaucld

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(291,2,im_index)),                  &
     &        meantaucld,                                               &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,291,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 291)"
            goto 9999
         End if

      End if

      If (sf(292,2)) Then

!   meanptop

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(292,2,im_index)),                  &
     &        meanptop,                                                 &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,292,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 292)"
            goto 9999
         End if

      End if

      If (sf(293,2)) Then

!   totalcldarea

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(293,2,im_index)),                  &
     &        totalcldarea,                                             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,293,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 293)"
            goto 9999
         End if

      End if

!   Cloud water mixing ratios
      item = 308  ! LS cloud liquid water mixing ratio
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        ls_qcl_rad,                                               &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 308)"//cmessage
         End if
      Endif
      item = 309  ! LS cloud ice water mixing ratio
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        ls_qcf_rad,                                               &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 309)"//cmessage
         End if
      Endif
      item = 310  ! Convective cloud liquid water mixing ratio
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cc_qcl_rad,                                               &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 310)"//cmessage
         End if
      Endif
      item = 311  ! Convective cloud ice water mixing ratio
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cc_qcf_rad,                                               &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 311)"//cmessage
         End if
      Endif
!   Cloud amounts
      item = 312  ! LS cloud fraction of grdbox seen by radiation.
                  ! Liquid
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        ls_cl_rad,                                                &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 312)"//cmessage
         End if
      Endif
      item = 313  ! LS cloud fraction of grdbox seen by radiation.
                  ! Ice
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        ls_cf_rad,                                                &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 313)"//cmessage
         End if
      Endif
      item = 314  ! CONV cloud fraction of grdbox seen by radiation.
                  ! Liquid
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cc_cl_rad,                                                &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 314)"//cmessage
         End if
      Endif
      item = 315  ! CONV cloud fraction of grdbox seen by radiation.
                  ! Ice
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cc_cf_rad,                                                &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 315)"//cmessage
         End if
      Endif

! Copy ISCCP diagnostics by looping over 7 levels in call to copydiag.
! This is because copydiag_3d cannot handle ISCCP levels.

      If (sf(270,2)) Then

!   isccp_cf

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(270,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf(1,1,k),                                          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,270,                                           &
     &        icode,cmessage)
        Enddo


         If (icode  >   0) then
            cmessage=": error in copydiag( item 270)"//cmessage
            goto 9999
         End if

      End if

      If (sf(271,2)) Then

!   isccp_cf_tau_0_to_p3

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(271,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_0_to_p3(1,1,k),                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,271,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 271)"//cmessage
            goto 9999
         End if

      End if

      If (sf(272,2)) Then

!   isccp_cf_tau_p3_to_1p3

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(272,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_p3_to_1p3(1,1,k),                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,272,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 272)"//cmessage
            goto 9999
         End if

      End if

      If (sf(273,2)) Then

!   isccp_cf_tau_1p3_to_3p6

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(273,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_1p3_to_3p6(1,1,k),                           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,273,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 273)"//cmessage
            goto 9999
         End if

      End if

      If (sf(274,2)) Then

!   isccp_cf_tau_3p6_to_9p4

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(274,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_3p6_to_9p4(1,1,k),                           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,274,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 274)"//cmessage
            goto 9999
         End if

      End if

      If (sf(275,2)) Then

!   isccp_cf_tau_9p4_to_23

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(275,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_9p4_to_23(1,1,k),                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,275,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 275)"//cmessage
            goto 9999
         End if

      End if

      If (sf(276,2)) Then

!   isccp_cf_tau_23_to_60

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(276,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_23_to_60(1,1,k),                             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,276,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 276)"//cmessage
            goto 9999
         End if

      End if

      If (sf(277,2)) Then

!   isccp_cf_tau_ge_60

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(277,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_ge_60(1,1,k),                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,277,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 277)"//cmessage
            goto 9999
         End if

      End if

!   Aerosol optical depth diagnostics
!   (loop on wavelength)

      If (sf(284,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(284,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_sulphate(1,1,k),                                      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,284,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 284)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(285,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(285,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_dust(1,1,k),                                          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,285,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 285)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(286,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(286,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_seasalt(1,1,k),                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,286,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 286)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(287,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(287,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_soot(1,1,k),                                          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,287,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 287)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(288,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(288,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_biomass(1,1,k),                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,288,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 288)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(289,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(289,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_biogenic(1,1,k),                                      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,289,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 289)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(295,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(295,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_ocff(1,1,k),                                          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,295,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 295)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(296,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(296,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_delta(1,1,k),                                         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,296,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 296)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

 9999 continue  ! exit point on error
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
      END SUBROUTINE diagnostics_lw
