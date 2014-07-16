
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_sw

      Subroutine diagnostics_sw(                                        &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels, cloud_levels             &
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
     &,                      surfsw, itoasw, solar_out_toa              &
     &,                      solar_out_clear, surface_down_flux         &
     &,                      surf_down_clr, surf_up_clr                 &
     &,                      SWsea, flux_below_690nm_surf               &
     &,                      photosynth_act_rad, dirpar_flux            &
     &,                      fl_solid_below_690nm_surf                  &
     &,                      fl_sea_below_690nm_surf                    &
     &,                      orog_corr, sol_bearing                     &
     &,                      f_orog                                     &
     &,                      slope_aspect, slope_angle                  &
     &,                      sw_net_land,sw_net_sice                    &
     &,                      T_incr_diagnostic                          &
     &,                      clear_hr                                   &
     &,                      flux_direct, flux_diffuse                  &
     &,                      net_flux_trop, up_flux_trop                &
     &,                      re_strat, wgt_strat, lwp_strat             &
     &,                      re_conv, wgt_conv                          &
     &,                      ntot_diag, strat_lwc_diag                  &
     &,                      so4_ccn_diag, cond_samp_wgt                &
     &,                      weighted_re, sum_weight_re                 &
     &,                      weighted_warm_re, sum_weight_warm_re       &
     &,                      Nc_diag, Nc_weight                         &
     &,                      sea_salt_film, sea_salt_jet                &
     &,                      salt_dim1, salt_dim2, salt_dim3            &
     &,                      cloud_extinction                           &
     &,                      cloud_weight_extinction                    &
     &,                      ls_cloud_extinction                        &
     &,                      ls_cloud_weight_extinction                 &
     &,                      cnv_cloud_extinction                       &
     &,                      cnv_cloud_weight_extinction                &
     &,                                                                 &
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
! Version   Date     Comment
! ---      -----     -------
! 5.0  30/11/99 Original version in UM. J-C Thil.
! 5.1  09/12/99 Add error trapping. Rick Rawlins
! ----     -------     -------
! 5.1  27/03/00 Add new diagnostics. J-C Thil.
! 5.1  28/02/00 Provide explicit model increments as STASH output
!               diagnostics. R Rawlins
! 5.2  14/11/00 Provide code for reactivated diagnostics.
!                                 (J. M. Edwards)
! 5.2  15/11/00 Add sea-salt aerosol diagnostics. A. Jones
! 5.3  25/04/01  Add diagnostics for coastal tiling.
!                                                   N. Gedney
! 5.4  03/08/01 Add PC2 cloud scheme diagnostics. D.R. Wilson
! 5.4  29/05/02 Add column-integrated cloud droplet and warm-cloud-only
!               satellite-view rE diagnostics.
!                                                 A. Jones
! 5.5  09/01/02 Add Extinction diagnostics. A.B.Keen/K.D.Williams
! 5.5  15/05/03 Correct cmessage for diagnostic 248.  A. Jones
! 6.2  02/03/06 Added diagnostics for total and direct component
!               of surface PAR flux.  M.G. Sanderson
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
     &, number_format                                                   &
                         ! switch controlling number format diagnostics
                         ! are written out in. See PP_WRITE for details.
     &, model_domain                                                    &
                         ! indicator as to model type, ie global, lam
     &, salt_dim1                                                       &
                         !
     &, salt_dim2                                                       &
                         ! Dimensions for sea-salt aerosol diagnostics.
     &, salt_dim3        !

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

      Real                                                              &
     &  timestep

! Primary Arrays used in all models

      Real                                                              &
     &  p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, T_n(row_length, rows, model_levels)                             &
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

     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                       1-off_y:rows+off_y, model_levels)          &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                   1-off_y:rows+off_y, model_levels)              &
     &, q(row_length,rows, wet_model_levels)                            &
     &, qcl(row_length, rows, wet_model_levels)                         &
     &, qcf(row_length, rows, wet_model_levels)

      Real                                                              &
     &  itoasw (row_length, rows)                                       &
     &, surfsw (row_length, rows)                                       &
     &, solar_out_toa (row_length, rows)                                &
     &, solar_out_clear (row_length, rows)                              &
     &, surface_down_flux (row_length, rows)                            &
     &, surf_down_clr (row_length, rows)                                &
     &, surf_up_clr (row_length, rows)                                  &
     &, SWsea(row_length, rows)                                         &
                                  ! Net short-wave absorbed by planet
     &, flux_below_690nm_surf(row_length, rows)                         &
     &, photosynth_act_rad(row_length, rows)                            &
     &, dirpar_flux(row_length, rows)                                   &
     &, fl_solid_below_690nm_surf(row_length, rows)                     &
     &, fl_sea_below_690nm_surf(row_length, rows)                       &
     &, sw_net_land(row_length, rows)                                   &
                                       !SW net local flux over land
     &, sw_net_sice(row_length, rows)                                   &
                                       !SW net local flux over sea-ice
     &, T_incr_diagnostic(row_length,rows,model_levels)                 &
     &, clear_hr(row_length, rows, model_levels)                        &
     &, flux_direct(row_length, rows, model_levels+1)                   &
     &, flux_diffuse(row_length, rows, model_levels+1)                  &
     &, net_flux_trop(row_length, rows)                                 &
     &, up_flux_trop(row_length, rows)                                  &
     &, re_strat(row_length, rows, cloud_levels)                        &
     &, wgt_strat(row_length, rows, cloud_levels)                       &
     &, lwp_strat(row_length, rows, cloud_levels)                       &
     &, re_conv(row_length, rows, cloud_levels)                         &
     &, wgt_conv(row_length, rows, cloud_levels)                        &
     &, ntot_diag(row_length, rows, cloud_levels)                       &
     &, strat_lwc_diag(row_length, rows, cloud_levels)                  &
     &, so4_ccn_diag(row_length, rows, cloud_levels)                    &
     &, cond_samp_wgt(row_length, rows, cloud_levels)                   &
     &, weighted_re(row_length, rows)                                   &
     &, sum_weight_re(row_length, rows)                                 &
     &, weighted_warm_re(row_length, rows)                              &
     &, sum_weight_warm_re(row_length, rows)                            &
     &, Nc_diag(row_length, rows)                                       &
     &, Nc_weight(row_length, rows)                                     &
     &, sea_salt_film(salt_dim1, salt_dim2, salt_dim3)                  &
     &, sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)                   &
     &, cloud_extinction(row_length, rows, cloud_levels)                &
     &, cloud_weight_extinction(row_length, rows, cloud_levels)         &
     &, ls_cloud_extinction(row_length, rows, cloud_levels)             &
     &, ls_cloud_weight_extinction(row_length, rows, cloud_levels)      &
     &, cnv_cloud_extinction(row_length, rows, cloud_levels)            &
     &, cnv_cloud_weight_extinction(row_length, rows, cloud_levels)

! Orography variables

      Real orog_corr(row_length, rows)   ! Orography correction factor
      Real sol_bearing(row_length,rows)  ! Local solar bearing
      Real f_orog(row_length,rows)       ! Extra SW surf flux
      Real slope_aspect(row_length,rows) ! Gridbox mean slope aspect
      Real slope_angle(row_length,rows)  ! Gridbox mean slope angle

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
     &  T_plus_T_inc(row_length, rows, model_levels)                    &
     &, heating_rate(row_length,rows,model_levels)                      &
     &, work_3d(row_length,rows,model_levels)

      Integer                                                           &
     &  i, j, k                                                         &
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      Character*80  cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_sw')

      Integer                                                           &
     &  im_index                                                        &
                        ! internal model index
     & ,item                                                            &
                        ! STASH item
     & ,sect            ! STASH section
      Parameter( sect = 1 ) ! for sw radiation

! External routines
      External                                                          &
     &  copydiag, copydiag_3d                                           &
     &  ,Ereport

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------

      If (sf(004,1)) Then

!L   T+T_inc

         Do k = 1, model_levels
            Do j = 1, rows
               Do i = 1, row_length
                  T_plus_T_inc(i,j,k) = T_n(i,j,k) + T_inc(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(004,1,im_index)),                &
     &        T_plus_T_inc,                                             &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,004,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,004,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage="Error in copydiag_3d( item 004)"
            goto 9999
         End if

      End if

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


      If (sf(201,1)) Then

!L   surfsw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(201,1,im_index)),surfsw,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,201,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if



      If (sf(207,1)) Then

!L   itoasw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(207,1,im_index)),itoasw,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,207,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 207)"
            goto 9999
         End if

      End if



      If (sf(208,1)) Then

!L   solar_out_toa

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(208,1,im_index)),solar_out_toa,    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,208,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
         End if

      End if


      If (sf(209,1)) Then

!L   solar_out_clear

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(209,1,im_index)),solar_out_clear,  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,209,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 209)"
            goto 9999
         End if

      End if



      If (sf(210,1)) Then

!L   surf_down_clr

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(210,1,im_index)),surf_down_clr,    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,210,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 210)"
            goto 9999
         End if

      End if



      If (sf(211,1)) Then

!L   surf_up_clr

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(211,1,im_index)),surf_up_clr,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,211,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 211)"
            goto 9999
         End if

      End if

      If (sf(235,1)) Then

!L   surface_down_flux

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(235,1,im_index)),surface_down_flux,&
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,235,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 235)"
            goto 9999
         End if

      End if


      If (sf(203,1)) Then

!L   SWsea

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(203,1,im_index)),SWsea,            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,203,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 203)"
            goto 9999
         End if

      End if


      If (sf(204,1)) Then

!L flux_below_690nm_surf

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(204,1,im_index)),                   &
     &        flux_below_690nm_surf,                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,204,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 204)"
            goto 9999
         End if

      End if
!
      If (sf(257,1)) Then

!L sw_net_land

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(257,1,im_index)),                   &
     &        sw_net_land,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,257,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 257)"
            goto 9999
         End if

      End if
!
      If (sf(258,1)) Then

!L sw_net_sice

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(258,1,im_index)),                   &
     &        sw_net_sice,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,258,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 258)"
            goto 9999
         End if

      End if
!
      If (sf(259,1)) Then

!L fl_solid_below_690nm_surf

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(259,1,im_index)),                   &
     &        fl_solid_below_690nm_surf,                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,259,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 259)"
            goto 9999
         End if

      End if
!
      If (sf(260,1)) Then

!L fl_sea_below_690nm_surf

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(260,1,im_index)),                   &
     &        fl_sea_below_690nm_surf,                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,260,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 260)"
            goto 9999
         End if

      End if
!


      If (sf(221,1)) Then

!L re_strat

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(221,1,im_index)),                &
     &        re_strat,                                                 &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,221,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,221,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 221)"
            goto 9999
         End if

      End if



      If (sf(223,1)) Then

!L wgt_strat

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223,1,im_index)),                &
     &        wgt_strat,                                                &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,223,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,223,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 223)"
            goto 9999
         End if

      End if



      If (sf(224,1)) Then

!L lwp_strat

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(224,1,im_index)),                &
     &        lwp_strat,                                                &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,224,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,224,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 224)"
            goto 9999
         End if

      End if



      If (sf(225,1)) Then

!L re_conv

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(225,1,im_index)),                &
     &        re_conv,                                                  &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,225,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,225,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 225)"
            goto 9999
         End if

      End if



      If (sf(226,1)) Then

!L wgt_conv

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(226,1,im_index)),                &
     &        wgt_conv,                                                 &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,226,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,226,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 226)"
            goto 9999
         End if

      End if

      If (sf(230,1)) Then
      
! flux_direct

         Call copydiag_3d(STASHwork(si(230,1,im_index)),                &
     &        flux_direct,                                              &
     &        row_length,rows,model_levels+1,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,230,1,im_index)),len_stlist,           &  
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,230,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 230)"
            goto 9999
         End if

      End if

      If (sf(231,1)) Then

! flux_diffuse

         Call copydiag_3d(STASHwork(si(231,1,im_index)),                &
     &        flux_diffuse,                                             & 
     &        row_length,rows,model_levels+1,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,231,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,231,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 231)"
            goto 9999
         End if

      End if

!
! SW heating =
!            sw radiation temperature increment per timestep / timestep

      If (icode <= 0 .and. sf(232,1)) Then

        Do k=1,model_levels
         Do j=1,rows
          Do i=1,row_length
            heating_rate(i,j,k) =  T_incr_diagnostic(i,j,k) /           &
     &                                  timestep
          Enddo ! i
         Enddo ! j
        Enddo ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(232,1,im_index)),                &
     &        heating_rate,                                             &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,232,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,232,                                           &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 232)"//cmessage
        End if

      Endif

      If (sf(233,1)) Then

!L clear_hr

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(233,1,im_index)),                &
     &        clear_hr,                                                 &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,233,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,233,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 233)"
            goto 9999
         End if

      End if


      If (sf(237,1)) Then

!L   net_flux_trop

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(237,1,im_index)),net_flux_trop,    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,237,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 237)"
            goto 9999
         End if

      End if


      If (sf(238,1)) Then

!L   up_flux_trop

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(238,1,im_index)),up_flux_trop,     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,238,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 238)"
            goto 9999
         End if

      End if



      If (sf(241,1)) Then

!L ntot_diag

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(241,1,im_index)),                &
     &        ntot_diag,                                                &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,241,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,241,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 241)"
            goto 9999
         End if

      End if



      If (sf(242,1)) Then

!L strat_lwc_diag

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(242,1,im_index)),                &
     &        strat_lwc_diag,                                           &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,242,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,242,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 242)"
            goto 9999
         End if

      End if



      If (sf(243,1)) Then

!L so4_ccn_diag

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(243,1,im_index)),                &
     &        so4_ccn_diag,                                             &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,243,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,243,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 243)"
            goto 9999
         End if

      End if



      If (sf(244,1)) Then

!L cond_samp_wgt

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(244,1,im_index)),                &
     &        cond_samp_wgt,                                            &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,244,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,244,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 244)"
            goto 9999
         End if

      End if



      If (sf(245,1)) Then

!L weighted_re

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(245,1,im_index)),                  &
     &        weighted_re,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,245,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 245)"
            goto 9999
         End if

      End if



      If (sf(246,1)) Then

!L sum_weight_re

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(246,1,im_index)),                  &
     &        sum_weight_re,                                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,246,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 246)"
            goto 9999
         End if

      End if



      If (sf(254,1)) Then

!L weighted_warm_re

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(254,1,im_index)),                  &
     &        weighted_warm_re,                                         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,254,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 254)"
            goto 9999
         End if

      End if



      If (sf(255,1)) Then

!L sum_weight_warm_re

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(255,1,im_index)),                  &
     &        sum_weight_warm_re,                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,255,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 255)"
            goto 9999
         End if

      End if



      If (sf(247,1)) Then

!L sea_salt_film

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(247,1,im_index)),                &
     &        sea_salt_film,                                            &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,247,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,247,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 247)"
            goto 9999
         End if

      End if



      If (sf(248,1)) Then

!L sea_salt_jet

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(248,1,im_index)),                &
     &        sea_salt_jet,                                             &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,248,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,248,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 248)"
            goto 9999
         End if

      End if



      If (sf(280,1)) Then

!L Nc_diag

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(280,1,im_index)),                  &
     &        Nc_diag,                                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,280,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 280)"
            goto 9999
         End if

      End if



      If (sf(281,1)) Then

!L Nc_weight

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(281,1,im_index)),                  &
     &        Nc_weight,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,281,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 281)"
            goto 9999
         End if

      End if


      If (sf(262,1)) Then

!   cloud_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(262,1,im_index)),                &
     &        cloud_extinction,                                         &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,262,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,262,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 262)"
            goto 9999
         End if

      End if

      If (sf(263,1)) Then

!   cloud_weight_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(263,1,im_index)),                &
     &        cloud_weight_extinction,                                  &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,263,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,263,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 263)"
            goto 9999
         End if

      End if

      If (sf(264,1)) Then

!   ls_cloud_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(264,1,im_index)),                &
     &        ls_cloud_extinction,                                      &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,264,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,264,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 264)"
            goto 9999
         End if

      End if

      If (sf(265,1)) Then

!   ls_cloud_weight_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(265,1,im_index)),                &
     &        ls_cloud_weight_extinction,                               &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,265,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,265,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 265)"
            goto 9999
         End if

      End if

      If (sf(266,1)) Then

!   cnv_cloud_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(266,1,im_index)),                &
     &        cnv_cloud_extinction,                                     &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,266,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,266,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 266)"
            goto 9999
         End if

      End if

      If (sf(267,1)) Then

!   cnv_cloud_weight_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(267,1,im_index)),                &
     &        cnv_cloud_weight_extinction,                              &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,267,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,267,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 267)"
            goto 9999
         End if

      End if


      If (sf(290,1)) Then

!L photosynth_act_rad (total PAR flux at surface)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(290,1,im_index)),                   &
     &        photosynth_act_rad,                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,290,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 290)"
            goto 9999
         End if

      End if

      If (sf(291,1)) Then

!L flux_direct_par (direct component of PAR flux at surface)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(291,1,im_index)),                   &
     &        dirpar_flux,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,291,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 291)"
            goto 9999
         End if

      End if


      If (sf(292,1)) Then

         !  sol_bearing

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(292,1,im_index)),sol_bearing,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,292,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 292)"
            goto 9999
         End if

      End if


      If (sf(293,1)) Then

         !  slope_aspect

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(293,1,im_index)),slope_aspect,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,293,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 293)"
            goto 9999
         End if

      End if


      If (sf(294,1)) Then

         !  slope_angle

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(294,1,im_index)),slope_angle,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,294,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 294)"
            goto 9999
         End if

      End if


      If (sf(295,1)) Then

         !  orog_corr

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(295,1,im_index)),orog_corr,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,295,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 295)"
            goto 9999
         End if

      End if


      If (sf(296,1)) Then

         !  f_orog

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(296,1,im_index)),f_orog,            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,296,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 296)"
            goto 9999
         End if

      End if


 9999 continue  ! exit point on error
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
      END SUBROUTINE diagnostics_sw

