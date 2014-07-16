
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Diagnostic output area for Large Scale Cloud (Section 9) Diagnostics.
! Subroutine Interface:
      Subroutine diagnostics_lscld(                                     &
     &                       row_length, rows, model_levels             &
     &,                      rhc_row_length, rhc_rows                   &
     &,                      wet_model_levels, boundary_layer_levels    &
     &,                      cloud_levels                               &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity, p_theta_levels               &
     &,                      p, r_theta_levels, r_rho_levels            &
     &,                      T_incr_diagnostic, q_incr_diagnostic       &
     &,                      qcl_incr_diagnostic                        &
     &,                      T, q, qcl, qcf                             &
     &,                      area_cloud_fraction, bulk_cloud_fraction   &
     &,                      p_star, rhcpt                              &
     &,                      combined_cloud, cca, ccb, cct              &
     &,                      n_cca_levels, L_murk, Aerosol, RHcrit      &
     &,                      l_mixing_ratio                             &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork                                                        &
     &  )
!
Use ac_diagnostics_mod, Only :  &
    cf_lsc, qcl_lsc

      Implicit None
!
! Purpose:
!          Calculates diagnostics associated with large scale cloud
! (section 9) for output through STASH system. This routine is
! positioned at the end of moist processes calculation, ie the end
! of the moist timestep, and is therefore the appropriate location
! for output of moisture related diagnostics.
!
! Method:
!  Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the calling
! routine. After minor calculations, each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing.
!
!  Diagnostics currently available: (in order calculated)
!
! STASH item (all section 9 )
!   4 Temperature on model levels after large scale cloud
! 181 T   increment over boundary layer + large scale cloud routines
! 182 q   increment over boundary layer + large scale cloud routines
! 183 qcl increment over boundary layer + large scale cloud routines
! 201 bulk cloud fraction
! 202 very low cloud amount
! 203 low      cloud amount
! 204 medium   cloud amount
! 205 high     cloud amount
! 206 qcl after large scale cloud. Note this is only needed to provide
!     assimilation increments for latent heat nudging. Otherwise should
!     use prognostic (0,254).
! 208 cloud base for cover >  0.1 octa kft
! 209 cloud base for cover >  1.5 octa kft
! 210 cloud base for cover >  2.5 octa kft
! 211 cloud base for cover >  3.5 octa kft
! 212 cloud base for cover >  4.5 octa kft
! 213 cloud base for cover >  5.5 octa kft
! 214 cloud base for cover >  6.5 octa kft
! 215 cloud base for cover >  7.9 octa kft
! 216 total cloud ; random overlap
! 217 total cloud ; max/random overlap
! 218 cloud fraction below 1000 ft asl
! 219 low cloud base ft asl
! 220 low cloud top ft asl
! 221 wet bulb freezing level
! 222 wet bulb temperature
! 223 total cloud top height kft
! 229 relative humidity on model levels. Note: percentage units.
! 226 layer cloud frequency
! 228 critical relative humidity on model levels
! 230 Visibility on model levels. Note that if murk aerosol is included
!     then vis depends on aerosol as well as humidity. Visibility is
!     available over all model levels but aerosol is only actively
!     mixed over boundary levels.
! 231 combined cloud amount on model levels
!
! Current Owner of Code: Owner of Section 9 (Large-Scale Cloud).
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
!  Global Variables:----------------------------------------------------
!     None.
!
! Arguments with Intent IN. ie: Input variables.

      LOGICAL ,INTENT(IN) ::                                            &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      INTEGER ,INTENT(IN) ::                                            &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, rhc_row_length                                                  &
                         ! row_length if diagnostic RHcrit ON, else 1.
     &, rhc_rows                                                        &
                         ! rows       if diagnostic RHcrit ON, else 1.
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, cloud_levels                                                    &
                         ! number of cloud model levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, boundary_layer_levels                                           &
                              ! number of boundary layer levels
     &, n_cca_levels     ! Number of levels for conv cloud
                         ! amount: 1 for 2D, nlevs for 3D.

      INTEGER ,INTENT(IN) ::                                            &
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

      LOGICAL ,INTENT(IN) ::                                            &
     &  L_murk                                                          &
                            ! murk aerosol included if T
     &, L_mixing_ratio      ! Use mixing ratio formulation

! Primary Arrays used in all models
      Real ,INTENT(IN) ::                                               &
     &  p_theta_levels(row_length, rows, model_levels)                  &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, T(row_length, rows, model_levels)                               &
     &, q(row_length,rows, wet_model_levels)                            &
     &, qcl(row_length, rows, wet_model_levels)                         &
     &, qcf(row_length, rows, wet_model_levels)                         &
     &, area_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, bulk_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, p_star(row_length, rows)                                        &
     &, RHCPT(rhc_row_length, rhc_rows, wet_model_levels)               &
     &, cca(row_length, rows, n_cca_levels)                             &
     &, T_incr_diagnostic(row_length, rows, model_levels)               &
     &, q_incr_diagnostic(row_length,rows, wet_model_levels)            &
     &, qcl_incr_diagnostic(row_length, rows, wet_model_levels)         &
     &, rhcrit(wet_model_levels)                                        &
                                    ! critical RH for cloud formation
     &, Aerosol(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
                                    ! murk aerosol
     &, combined_cloud(row_length, rows, wet_model_levels)
!                                        Mixed CCA and CF per gridbox.
!
      INTEGER  ,INTENT(IN) ::                                           &
     &  ccb(row_length, rows)                                           &
     &, cct(row_length, rows)

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
!         cconsts defines LOW_BOT_LEVEL to HIGH_TOP_LEVEL cloud splits.
! CLSCHM3A defines reference numbers for cloud schemes in two-stream
! radiation code.

      ! maximum/random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_MAX=2

      ! random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_RANDOM=4

      ! maximum overlap in a column model
      INTEGER,PARAMETER:: IP_CLOUD_COLUMN_MAX=3

      ! clear column
      INTEGER,PARAMETER:: IP_CLOUD_CLEAR=5

      ! mixed column with split between  convective and layer cloud.
      INTEGER,PARAMETER:: IP_CLOUD_TRIPLE=6

      ! Coupled overlap with partial correlation of cloud
      INTEGER,Parameter:: IP_cloud_part_corr=7

      ! Coupled overlap with partial correlation of cloud
      ! with a separate treatment of convective cloud
      INTEGER,Parameter:: IP_cloud_part_corr_cnv=8

! CLSCHM3A end
!         clschm3a defines reference numbers for cloud overlap choices.
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
! C_KT_FT start

      REAL,PARAMETER:: KT2MS=1852.0/3600.0 ! Knots to m/s conversion
      REAL,PARAMETER:: FT2M =0.3048        ! Feet to meters conversion

! C_KT_FT end
! C_LOWCLD contains declarations for ceiling and threshold
! constants for lowest cloud layer diagnostics.
      ! Max height asl for 'low' cloud (1000ft)
      REAL,PARAMETER:: STR_CEIL=1000.0

      ! Cloud fraction threshold for low cloud.
      REAL,PARAMETER:: CLOUD_THRESHOLD=0.05
! C_LOWCLD end
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
!*L------------------COMDECK C_A----------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Replace variable A by more meaningful name for
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Convert to Fixed/Free format. P. Selwood

      ! Mean radius of Earth in m.
      Real, Parameter  :: Earth_Radius = 6371000.

!*----------------------------------------------------------------------
!
! arguments with intent in/out. ie: input variables changed on output.
!
! Diagnostics info
      Real                                                              &
     & STASHwork(*)     ! STASH workspace
!
!  Local parameters and other physical constants------------------------
!
      Integer sect ! STASH section for diagnostics
      Parameter ( sect = 9 )
!      ( LS Cloud diagnostics are output under STASH Section 09 )
!
      Character(*) RoutineName
      Parameter  ( RoutineName='diagnostics_lscld')
!
!  Local scalars--------------------------------------------------------
!
      Integer                                                           &
     &  i, j, k, l, ji                                                  &
     &, kinvert                                                         &
                                ! vertical index for inverted arrays.
     &, item                                                            &
                                ! STASH item of individual diagnostics.
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      REAL, PARAMETER:: vis_probability=0.5 ! median visibility diag

      Character*80  cmessage
!
      Integer                                                           &
     &  im_index                                                        &
                        ! internal model index
!                  ( LS Cloud is part of the Atmosphere internal model )
     &, nclds      ! Number of radiation cloud levels ( <=  wet levels)
!  Local scalars required to calculate Diagnostics 09208 - 09215
      INTEGER, PARAMETER :: noktas=8
      REAL,    PARAMETER :: roktas=8.0
      REAL    :: cloud_base_height         ! Height of cloud base (m)
      REAL    :: combined_cloud_in_oktas   ! Combined cloud in oktas

!  Local dynamic arrays required to calculate Diagnostics 09208 - 09215
      INTEGER ::                                                        &
     & cloud_base_level(row_length, rows, noktas)
      REAL ::                                                           &
     & cloud_base(row_length, rows, noktas)
      REAL :: oktas(noktas)
      DATA oktas/0.1,1.5,2.5,3.5,4.5,5.5,6.5,7.9/

!  Local scalars required to calculate Diagnostics 09218, 09219, 09220
      REAL ::                                                           &
     & pu,pl                                                            &
                                       ! Upper and lower half
                                       !       level pressure
     &,pt                                                               &
                                       ! Pressure thickness accumulator
     &,ft                                                               &
                                       ! Cloud fraction accumulator
     &,dp                                                               &
                                       ! Later pressure thickness
     &,h_asl                                                            &
                                       ! Layer base height asl
     &,h_asln                                                           &
                                       ! Layer top height asl
     &,fr                              ! Layer fraction below ceiling

      REAL, PARAMETER ::                                                &
     & str_ceilm = str_ceil * ft2m     ! STR_CEIL in metres

!  Local scalars required to calculate Diagnostic 09221 and 09222
      LOGICAL, PARAMETER ::                                             &
     & l_potential=.false.             ! Wet bulb temperature required
                                       ! from subroutine ThetaW
      REAL ::                                                           &
     & frac

!  Local scalars required to calculate Diagnostic 09223
      INTEGER ::                                                        &
     & cld_top_lev                     ! Top layer in cloud
      REAL, PARAMETER::                                                 &
     & thresh = 0.0627,                                                 &
     & m_to_kft = (1./ft2m)*0.001
!
!  Local dynamic arrays-------------------------------------------------
       Real                                                             &
     &  work2d_1(row_length, rows)                                      &
                                       ! Single-level work array (cloud)
     &, work2d_2(row_length, rows)                                      &
     &, work2d_3(row_length, rows)                                      &
     &, work3d_1(row_length, rows, wet_model_levels)
!                                        Full work array (cloud/moist).
!
!  External subroutine calls: ------------------------------------------
      External                                                          &
     &  R2_calc_total_cloud_cover                                       &
     &, QSAT                                                            &
     &, copydiag, copydiag_3d                                           &
     &, Ereport                                                         &
     &, Thetaw
!- End of Header
!
! ==Main Block==--------------------------------------------------------
      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)
!
! ----------------------------------------------------------------------
! DIAG.09004 Copy T from Main Cloud to stashwork
! ----------------------------------------------------------------------
      item = 4
! Diag09004_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),T,         &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
         If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(T Main Cloud)"
         End if
      End if  ! Diag09004_if1
!
! ----------------------------------------------------------------------
! DIAG.09010 Copy q to stashwork
! ----------------------------------------------------------------------
      If (.false.) Then
! Prevent access to diagnostic which duplicates 00010 and thus redundant
! Pending complete removal of code at UM Version 5.3.
      item = 10
! Diag09010_if1:
!     If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),q,         &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
!#include <argppx/argppx.h>
     &      icode,cmessage)

         If (icode >  0) Then
            cmessage="cld_ctl  : error in copydiag_3d(q)"
         End if
      End if  ! Diag09010_if1
!
! ----------------------------------------------------------------------
! DIAG.09181 Copy T INC: bdy layer + ls cld to stashwork
! ----------------------------------------------------------------------
      item = 181
! Diag09181_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       T_incr_diagnostic,                                         &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(T INC: bl + ls cld)"
        End if
      End if  ! Diag09181_if1
!
! ----------------------------------------------------------------------
! DIAG.09182 Copy q INC: bdy layer + ls cld to stashwork
! ----------------------------------------------------------------------
      item = 182
! Diag09182_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       q_incr_diagnostic,                                         &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Q INC: bl + ls cld)"
        End if
      End if  ! Diag09182_if1
!
! ----------------------------------------------------------------------
! DIAG.09183 Copy qCL INC: bdy layer + ls cld to stashwork
! ----------------------------------------------------------------------
      item = 183
! Diag09183_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
        Call copydiag_3d ( stashwork(si(item,sect,im_index)),           &
     &       qcl_incr_diagnostic,                                       &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
        cmessage="cld_ctl  : error in copydiag_3d(Qcl INC: bl + ls cld)"
        End if
      End if  ! Diag09183_if1
!
! ----------------------------------------------------------------------
! DIAG.09201 Copy Bulk cloud fraction to stashwork
! ----------------------------------------------------------------------
! Needed for Assimilation, otherwise duplicates 00266 output (preferred)
      item = 201
! Diag09201_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then

      If (.not.allocated(cf_lsc)) then
        Allocate ( cf_lsc(row_length*rows,wet_model_levels) )
      End If
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        bulk_cloud_fraction,                                      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
! Copy bulk_cloud_fraction into cf_lsc for ac_diagnostics_mod

      do k = 1,wet_model_levels
        do j = 1,rows
          do i = 1,row_length
            ji = (j-1)*row_length+i
            cf_lsc(ji,k) = bulk_cloud_fraction(i,j,k)
          end do
        end do
      end do

        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Bulk CF Main Cloud)"
        End if
      End if  ! Diag09201_if1
!
! ----------------------------------------------------------------------
! DIAG.09202 Copy Very Low Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 202
! Diag09202_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09202_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i, j) = area_cloud_fraction(i, j, 1)
          End Do
        End Do  ! Diag09203_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09202_do2:
        Do k = 1, MAX((LOW_BOT_LEVEL - 1), 1)
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
     &             MAX(work2d_1(i, j), area_cloud_fraction(i, j, k))
            End Do
          End Do
        End Do  ! Diag09202_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(LOW Cloud Amount)"
         End if
      End if  ! Diag09203_if1
!
! ----------------------------------------------------------------------
! DIAG.09203 Copy Low Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 203
! Diag09203_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09203_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i, j) = area_cloud_fraction(i, j, LOW_BOT_LEVEL)
          End Do
        End Do  ! Diag09203_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09203_do2:
        Do k = LOW_BOT_LEVEL + 1, LOW_TOP_LEVEL
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
     &             MAX(work2d_1(i, j), area_cloud_fraction(i, j, k))
            End Do
          End Do
        End Do  ! Diag09203_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(LOW Cloud Amount)"
         End if
      End if  ! Diag09203_if1
!
! ----------------------------------------------------------------------
! DIAG.09204 Copy Medium Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 204
! Diag09204_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09204_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i, j) = area_cloud_fraction(i, j, MED_BOT_LEVEL)
          End Do
        End Do  ! Diag09204_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09204_do2:
        Do k = MED_BOT_LEVEL + 1, MED_TOP_LEVEL
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
     &             MAX(work2d_1(i, j), area_cloud_fraction(i, j, k))
            End Do
          End Do
        End Do  ! Diag09204_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(MEDIUM Cloud Amount)"
        End if
      End if  ! Diag09204_if1
!
! ----------------------------------------------------------------------
! DIAG.09205 Copy High Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 205
! Diag09205_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09205_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i, j) = area_cloud_fraction(i, j, HIGH_BOT_LEVEL)
          End Do
        End Do  ! Diag09205_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09205_do2:
        Do k = HIGH_BOT_LEVEL + 1, HIGH_TOP_LEVEL
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
     &             MAX(work2d_1(i, j), area_cloud_fraction(i, j, k))
            End Do
          End Do
        End Do  ! Diag09205_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(HIGH Cloud Amount)"
        End if
      End if  ! Diag09205_if1
!
! ----------------------------------------------------------------------
! DIAG.09206 Copy cloud liquid condensate to stashwork
! ----------------------------------------------------------------------
! Needed for Assimilation, otherwise duplicates 00254 output (preferred)
      item = 206
! Diag09206_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
      If (.not.allocated(qcl_lsc)) Then
        Allocate ( qcl_lsc(row_length*rows,wet_model_levels) )
      End If

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(stashwork(si(item,sect,im_index)),qcl,        &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

! Copy qcl to qcl_lsc for ac_diagnostics_mod

      do k = 1,wet_model_levels
        do j = 1,rows
          do i = 1,row_length
            ji = (j-1)*row_length+i
            qcl_lsc(ji,k) = qcl(i,j,k)
          end do
        end do
      end do

         If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Qcl Main Cloud)"
         End if
      End if  ! Diag09206_if1
!
! ----------------------------------------------------------------------
! DIAG.09207 Copy cloud frozen condensate to stashwork
! ----------------------------------------------------------------------
      If (.false.) Then
! Prevent access to diagnostic after repositioning of routine until a
! final decision is made on if/how to ingest it. Timetable UM5.1 -> 5.2
      item = 207
! Diag09207_if1:
!     If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(stashwork(si(item,sect,im_index)),qcf,        &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage="cld_ctl  : error in copydiag_3d(cloud ice)"
         End if
      End if  ! Diag09207_if1
!
! ----------------------------------------------------------------------
! Find cloud base for pre defined cloud cover threshholds
! for Diagnostics 09208 - 09215.
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. ( sf(208,sect) .OR.                       &
     &                          sf(209,sect) .OR.                       &
     &                          sf(210,sect) .OR.                       &
     &                          sf(211,sect) .OR.                       &
     &                          sf(212,sect) .OR.                       &
     &                          sf(213,sect) .OR.                       &
     &                          sf(214,sect) .OR.                       &
     &                          sf(215,sect))) THEN
! Initialise output arrays
        cloud_base(:,:,:)=RMDI
        cloud_base_level(:,:,:)=IMDI

! Set cloud base to model levels
        DO l=1,noktas                 ! Loop over threshholds
          DO k=1,wet_model_levels     ! Loop over levels
            DO j=1,rows               ! Loop over rows
              DO i=1,row_length       ! Loop over points
                kinvert = wet_model_levels + 1 - k
  ! Convert to oktas
                combined_cloud_in_oktas=combined_cloud(i,j,kinvert)     &
     &                                  *roktas
  ! Calculate cloud_base_level
                IF(combined_cloud_in_oktas >  oktas(l))THEN
                  IF(cloud_base_level(i,j,l) <  0)THEN
                    cloud_base_level(i,j,l) = k
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO

! Convert level numbers to heights (M converted to Kft)
        DO l=1,noktas                 ! Loop over threshholds
          DO j=1,rows                 ! Loop over rows
            DO i=1,row_length         ! Loop over points
              IF(cloud_base_level(i,j,l) > 0)THEN
                cloud_base_height =                                     &
     &                 r_rho_levels(i,j,cloud_base_level(i,j,l))        &
     &                 - earth_radius
                cloud_base(i,j,l)=cloud_base_height*m_to_kft
              END IF
            END DO
          END DO
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09208 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 0.1 okta
! ----------------------------------------------------------------------
      item = 208
! Diag09208_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,1),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 0.1 okta)"
        END IF
      END IF  ! Diag09208_if1
!
! ----------------------------------------------------------------------
! DIAG.09209 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 1.5 okta
! ----------------------------------------------------------------------
      item = 209
! Diag09209_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,2),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 1.5 okta)"
        END IF
      END IF  ! Diag09209_if1
!
! ----------------------------------------------------------------------
! DIAG.09210 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 2.5 okta
! ----------------------------------------------------------------------
      item = 210
! Diag09210_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,3),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 2.5 okta)"
        END IF
      END IF  ! Diag09210_if1
!
! ----------------------------------------------------------------------
! DIAG.09211 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 3.5 okta
! ----------------------------------------------------------------------
      item = 211
! Diag09211_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,4),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 3.5 okta)"
        END IF
      END IF  ! Diag09211_if1
!
! ----------------------------------------------------------------------
! DIAG.09212 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 4.5 okta
! ----------------------------------------------------------------------
      item = 212
! Diag09212_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,5),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 4.5 okta)"
        END IF
      END IF  ! Diag09212_if1
!
! ----------------------------------------------------------------------
! DIAG.09213 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 5.5 okta
! ----------------------------------------------------------------------
      item = 213
! Diag09213_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,6),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 5.5 okta)"
        END IF
      END IF  ! Diag09213_if1
!
! ----------------------------------------------------------------------
! DIAG.09214 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 6.5 okta
! ----------------------------------------------------------------------
      item = 214
! Diag09214_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,7),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 6.5 okta)"
        END IF
      END IF  ! Diag09214_if1
!
! ----------------------------------------------------------------------
! DIAG.09215 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 7.9 okta
! ----------------------------------------------------------------------
      item = 215
! Diag09215_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,8),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 7.9 okta)"
        END IF
      END IF  ! Diag09215_if1
!
! ***** Code matches definition in COMMON RADIATION Section 70. ****
      nclds = MIN(cloud_levels, wet_model_levels)
!
! ----------------------------------------------------------------------
! DIAG.09216 Copy Total Cloud Amount RANDOM Overlap to stashwork
! ----------------------------------------------------------------------
      item = 216
! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
! Diag09216_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: r2_calc_total_cloud_cover
        Call R2_calc_total_cloud_cover(                                 &
     &         row_length*rows, wet_model_levels, nclds                 &
     &       , IP_CLOUD_MIX_RANDOM, combined_cloud(1,1,1), work2d_1     &
     &       , row_length*rows, wet_model_levels                        &
     &        )
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="cld_ctl  : error in copydiag(total cloud Random)"
         Endif
      Endif  ! Diag09216_if1
!
! ----------------------------------------------------------------------
! DIAG.09217 Copy Total Cloud Amount MAX/RANDOM Overlap to stashwork
! ----------------------------------------------------------------------
      item = 217
! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
! Diag09217_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: r2_calc_total_cloud_cover
        Call R2_calc_total_cloud_cover(                                 &
     &         row_length*rows, wet_model_levels, nclds                 &
     &       , IP_CLOUD_MIX_MAX, combined_cloud(1,1,1), work2d_1        &
     &       , row_length*rows, wet_model_levels                        &
     &        )
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(total cloud Max/Rand)"
        Endif
      Endif  ! Diag09217_if1
! ----------------------------------------------------------------------
! Find Cloud Fraction in air &lt; STR_CEIL, Low Cloud Base and Top
! for Diagnostics 09218 and/or 09219 and/or 09220
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. ( sf(218,sect) .OR.                       &
     &                          sf(219,sect) .OR.                       &
     &                          sf(220,sect))) THEN
        DO j=1,rows               ! Loop over rows
          DO i=1,row_length       ! Loop over points
            pt=0.0
            ft=0.0
            work2d_1(i,j)=RMDI    ! Initialise work arrays
            work2d_2(i,j)=RMDI
            work2d_3(i,j)=RMDI
            h_asln=r_theta_levels(i,j,0) - earth_radius
            pu=p_star(i,j)
            DO k=1,min(wet_model_levels,model_levels-1)
              h_asl=h_asln
              h_asln=r_rho_levels(i,j,k+1) - earth_radius

! Check if have not already found low cloud base
              IF(work2d_2(i,j) == RMDI) THEN
! IF not, and cloud cover in this layer &gt; threshhold
                IF(area_cloud_fraction(i,j,k) >= cloud_threshold) THEN
! THEN call the bottom of this layer the base of low cloud
                  work2d_2(i,j) = h_asl / ft2m
                END IF
              END IF

! Check if already found low cloud base but not top
              IF((work2d_2(i,j) /= RMDI).AND.(work2d_3(i,j) == RMDI))   &
     &          THEN
! IF not, and cloud cover in this layer &lt; threshhold
                IF(area_cloud_fraction(i,j,k) <= cloud_threshold) THEN
! THEN call the bottom of this layer the top of low cloud
                  work2d_3(i,j) = h_asl / ft2m
                END IF
              END IF

! IF bottom of layer is below low cloud ceiling (1000ft)
               IF(h_asl  <   str_ceilm) THEN
! Calculate top and bottom layer pressures
                 pl = pu
                 pu = p(i,j,k+1)
! And accumulate presure thickness and pressure weighted cloud amount
                 dp = pu - pl
! IF whole layer below ceiling, simply accumulate whole layer.
                 IF(h_asln  <   str_ceilm) THEN
                   pt = pt + dp
                   ft = ft + dp * area_cloud_fraction(i,j,k)
                 ELSE
! Otherwise height interpolate
                  fr = (str_ceilm - h_asl ) / (h_asln - h_asl )
                  pt = pt + dp * fr
                  ft = ft + dp * area_cloud_fraction(i,j,k) * fr
! And set results
                  work2d_1(i,j) = ft / pt

                 END IF
               END IF
             END DO               ! k over max
           END DO                 ! END loop over points
         END DO                   ! END loop over rows
       END IF
!
! ----------------------------------------------------------------------
! DIAG.09218 Copy Cloud fraction below 1000ft ASL to stashwork.
! ----------------------------------------------------------------------
      item = 218
! Diag09218_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_1,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     & "cld_ctl  : error in copydiag(Cloud fraction below 1000ft ASL)"
        END IF
      END IF  ! Diag09218_if1
!
! ----------------------------------------------------------------------
! DIAG.09219 Copy Low Cloud Base to stashwork.
! ----------------------------------------------------------------------
      item = 219
! Diag09219_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_2,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Low cloud base)"
        END IF
      END IF  ! Diag09219_if1
!
! ----------------------------------------------------------------------
! DIAG.09220 Copy Low Cloud Top to stashwork.
! ----------------------------------------------------------------------
      item = 220
! Diag09220_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_3,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Low cloud top)"
        END IF
      END IF  ! Diag09220_if1
!
! ----------------------------------------------------------------------
! Calculate wet bulb temperature for Diagnostics 09221 and/or 09222
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. ( sf(221,sect) .OR.                       &
     &                          sf(222,sect))) THEN

        DO k = 1, wet_model_levels

! Initialize 2d work arrays
          DO j = 1, rows
            DO i = 1, row_length
              work2d_1(i,j) = t(i,j,k)
              work2d_2(i,j) = q(i,j,k)
              work2d_3(i,j) = p_theta_levels(i,j,k)
            END DO
          END DO

! Calculate tw
! DEPENDS ON: thetaw
          Call ThetaW(row_length*rows,                                  &
     &                work2d_1, work2d_2, work2d_3, l_potential,        &
     &                work3d_1(1,1,k))
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09222 Copy Wet Bulb Temperature to stashwork.
! ----------------------------------------------------------------------
      item = 222
! Diag09222_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d(stashwork(si(item,sect,im_index)), work3d_1,   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag_3d(Wet Bulb Temperature)"
        END IF
      END IF  ! Diag09222_if1
!
! ----------------------------------------------------------------------
! Calculate wet bulb freezing level for Diagnostic 09221
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND.  sf(221,sect)) THEN
        DO j = 1, rows                 ! Loop over rows
          DO i = 1, row_length         ! Loop over points
            DO k = 1, wet_model_levels ! Loop over all wet-levels
              IF (work3d_1(i,j,k)  /=  rmdi) THEN
                IF (work3d_1(i,j,k)  <=  zerodegc) THEN
                  IF ( k  ==  1) THEN
                    work2d_1(i,j)=r_theta_levels(i,j,k) -               &
     &                            earth_radius
                  ELSE
                    frac = (zerodegc - work3d_1(i,j,k-1))/              &
     &                     (work3d_1(i,j,k)-work3d_1(i,j,k-1))
                    work2d_1(i,j)=r_theta_levels(i,j,k)*frac +          &
     &                            r_theta_levels(i,j,k-1)*(1.0-frac)    &
     &                            - earth_radius
                  END IF
                  EXIT
                END IF
              ELSE
                work2d_1(i,j) = rmdi
              END IF
            END DO
          END DO
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09221 Copy Wet Bulb Freezing Level to stashwork.
! ----------------------------------------------------------------------
      item = 221
! Diag09221_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_1,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Wet Bulb Freezing Level)"
        END IF
      END IF  ! Diag09221_if1

!
! ----------------------------------------------------------------------
! Calculate Total Cloud Top Height for Diagnostic 09223
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. sf(223,sect)) THEN
        DO j = 1, rows                   ! Loop over rows
          DO i = 1, row_length           ! Loop over points
            cld_top_lev = IMDI
            DO k = 1, wet_model_levels   ! Loop over all wet-levels
              IF (combined_cloud(i,j,k)  >   thresh) THEN
                cld_top_lev = wet_model_levels + 1 - k
              END IF
              IF (cld_top_lev  >   0) EXIT
            END DO
            IF (cld_top_lev  >   0) THEN
              work2d_1(i,j) = r_theta_levels(i,j,cld_top_lev)           &
     &                        - earth_radius
! Convert to Kiloft
              work2d_1(i,j) = work2d_1(i,j) * m_to_kft
            ELSE
              work2d_1(i,j) = RMDI
            END IF
          END DO
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09223 Copy Total Cloud Top Height to stashwork.
! ----------------------------------------------------------------------
      item = 223
! Diag09223_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_1,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Total cloud top height)"
        END IF
      END IF  ! Diag09223_if1
!
! ----------------------------------------------------------------------
! DIAG.09229 Copy Relative Humidity to stashwork
! ----------------------------------------------------------------------
      item = 229
! Diag09229_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Calculate Qsat then work out Relative Humidity.
! Diag09229_do1:
        Do k = 1, wet_model_levels
! DEPENDS ON: qsat_mix
          Call qsat_mix(work2d_1,T(1,1,K),p_theta_levels(1,1,k),        &
     &              row_length*rows,l_mixing_ratio)
!
          Do j = 1, rows
            Do i = 1, row_length
              ! relative humidity in per cent
              work3d_1(i, j, k) = q(i, j, k) / work2d_1(i, j) *100.0

!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
              If ( work3d_1(i, j, k) < 0.0) Then
                  work3d_1(i, j, k) = 0.0
              Endif

            End Do
          End Do
!
        End Do  ! Diag09229_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(RH Main Cloud)"
        End if
      End if  ! Diag09229_if1
!
! ----------------------------------------------------------------------
! DIAG.09226 Copy Layer Cloud Frequency to stashwork.
! ----------------------------------------------------------------------
      item = 226
! Diag09226_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Freq_k_do1:
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              If (area_cloud_fraction(i, j, k)  <=  0.)  Then
                work3d_1(i, j, k) = 0.
              Else
                work3d_1(i, j, k) = 1.
              End if
            End Do
          End Do
        End Do  ! Freq_k_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),work3d_1,   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Layer Cloud Freq)"
        End if
      End if  ! Diag09226_if1
!
! ----------------------------------------------------------------------
! DIAG.09228 Copy Critical Relative Humidity to stashwork.
! ----------------------------------------------------------------------
      item = 228
! Diag09228_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!     STASH will only permit this diagnostic to be chosen if 3D RHcrit
!     diagnosis is active. So RHCPT should be (row_length,rows,wet_ml).
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),RHCPT,      &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Critical RH)"
        End if
      End if  ! Diag09228_if1
!
! ----------------------------------------------------------------------
! DIAG.09230 Calculate visibility on model levels and copy to stashwork
! ----------------------------------------------------------------------
      item = 230
! Note that the following calculates the visibility for all model levels
! and then extracts just those levels that have been requested. If the
! diagnostic is normally called for a small subset of levels it will be
! more efficient to introduce a call to set_levels_list and calculate
! only those levels that are needed.
! Diag09230_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then

        do k=1,wet_model_levels

          ! visibility is calculated using prescribed critical RH
          ! (not diagnostic RH)
! DEPENDS ON: visbty
          Call Visbty(                                                  &
     &     p_theta_levels(1,1,k), T(1,1,k), q(1,1,k)                    &
                                                           !INPUT
     &     ,qcl(1,1,k), qcf(1,1,k)                                      &
                                                           !INPUT
     &     ,Aerosol(1:row_length,1:rows,k)                              &
                                                           !INPUT
     &     ,vis_probability, RHcrit(k), L_murk                          &
                                                           !INPUT
     &     ,row_length * rows                                           &
                                                           !INPUT
     &     ,work3d_1(1,1,k))                               !OUTPUT

        End do ! k   wet_model_levels

! DEPENDS ON: copydiag_3d
        Call copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)


      End if  ! Diag09230_if1

! ----------------------------------------------------------------------
! DIAG.09231 Copy Combined Cloud fraction to stashwork
! ----------------------------------------------------------------------
      item = 231
! Diag09231_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Diag09231_do1:
        Do k = 1, wet_model_levels
! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
!     Combined_cloud is calculated this way but re-inverted for STASH.
!
          kinvert = wet_model_levels+1-k
          Do j = 1, rows
            Do i = 1, row_length
              work3d_1(i, j, k) = combined_cloud(i,j,kinvert)
            End Do
          End Do
        End Do  ! Diag09231_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)), work3d_1,  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Combined Cld On Lv)"
        End if
      End if  ! Diag09231_if1
!
! ----------------------------------------------------------------------
      If (icode  /=  0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      End if
!
      Return
      END SUBROUTINE diagnostics_lscld
