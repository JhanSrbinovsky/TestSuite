
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine NI_filter_incs_Ctl

      Subroutine NI_filter_incs_Ctl(                                    &
     &                      theta, theta_star, R_u, R_v, R_w,           &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      r_theta_levels, r_rho_levels,               &
     &                      r_at_u, r_at_v,                             &
     &                      delta_lambda, delta_phi,                    &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      cos_v_latitude, sec_v_latitude,             &
     &                      north_lat_limit, south_lat_limit,           &
     &                      polar_filter_coefficient,                   &
     &                      polar_filter_n_sweeps, step_per_sweep,      &
     &                      polar_start_lat_limit,                      &
     &                      max_filter_rows, u_sweeps, v_sweeps,        &
     &                      global_u_filter, global_v_filter,           &
     &                      u_begin, u_end, v_begin, v_end,             &
     &                      diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
     &                      diff_coeff_thermo, diff_coeff_wind,         &
     &                      diff_order_thermo, diff_order_wind,         &
     &                      horizontal_level, me,                       &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      n_proc, n_procx, n_procy, l_datastart,      &
     &                      neighbour, at_extremity, model_domain,      &
     &                      proc_row_group, proc_col_group,             &
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
     &                      STASHwork)

! Purpose:
!          Call original polar filter code or new code
!          which also does horizontal diffusion as a filter
!
! Original Programmer: T. Davies
! Current code owner: T. Davies
!
! History:
! Version    Date      Comment
! ----     -------     -------
!  6.1     04/08/04    Code Introduced                     Terry Davies
!  6.2     25/12/05    Changes to filtering               Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

!include "parvars.h"
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

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  L_filter_incs                                                   &
                          ! switch for filtering of increments
     &, Ltimer                                                          &
     &, L_polar_filter_incs                                             &
                             ! switch for polar filtering of increments
     &, L_pofil_new                                                     &
                          ! new filter/diffusion switch
     &, L_pfcomb                                                        &
                          ! combined polar filter/diffusion active
     &, L_pftheta                                                       &
                          ! switch for polar filter for theta
     &, L_pfw                                                           &
                          ! switch for polar filter for w
     &, L_pfuv                                                          &
                          ! switch for polar filter for horizontal winds
     &, L_diff_incs       ! switch for horiz. diffusion of increments

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, l_datastart(3)                                                  &
                             ! First gridpoints held by this processor
     &, me                                                              &
                 !  this processor
     &, model_domain

      Real                                                              &
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
     &  cos_theta_longitude(row_length,rows)                            &
                                             ! cos of longitude
     &, sin_theta_longitude(row_length,rows)                            &
                                             ! sin of longitude
     &, sin_theta_latitude (row_length, rows)                           &
     &, sin_v_latitude (row_length, n_rows)                             &
     &, cos_theta_latitude(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y )                         &
     &, sec_theta_latitude(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y )                         &
     &, cos_v_latitude(1-off_x:row_length+off_x, 1-off_y:rows+off_y )   &
     &, sec_v_latitude(1-off_x:row_length+off_x, 1-off_y:rows+off_y )   &
!   diffusion coefficient for u/th rows
     &, diff_coeff_u(1-off_x:row_length+off_x,1-off_y:rows+off_y )      &
!   diffusion coefficient for v rows
     &, diff_coeff_v(1-off_x:row_length+off_x,1-off_y:n_rows+off_y )

      Real                                                              &
           ! limits for filtering (in radians)
     &  north_lat_limit                                                 &
                         ! latitude above which 1-2-1 filter is applied
     &, south_lat_limit                                                 &
                         ! latitude below which 1-2-1 filter is applied
     &, polar_filter_coefficient                                        &
                                     !   filter coefficient
     &, step_per_sweep                                                  &
                         ! amount in radians to increment start latitude
                         ! by per sweep
     &, polar_start_lat_limit                                           &
                              ! max latitude at which filter can start
     &, diff_coeff_wind                                                 &
                           ! NS diffusion coeff for winds
     &, diff_coeff_thermo                                               &
                           ! NS diffusion coeff for thermo
     &, diff_coeff_phi    ! NS diffusion coeff


      Integer                                                           &
     &  polar_filter_n_sweeps                                           &
                                     ! number of sweeps of filter to do
     &, max_filter_rows                                                 &
                           ! array size for u_begin etc
     &, global_u_filter                                                 &
                         ! sweep control for 1-2-1 filter
     &, global_v_filter                                                 &
                         ! sweep control for 1-2-1 filter
     &, u_begin(0:max_filter_rows)                                      &
                                    ! row pointers for 1-2-1 filter
     &, u_end(0:max_filter_rows)                                        &
                                    ! row pointers for 1-2-1 filter
     &, v_begin(0:max_filter_rows)                                      &
                                    ! row pointers for 1-2-1 filter
     &, v_end(0:max_filter_rows)                                        &
                                    ! row pointers for 1-2-1 filter
     &, u_sweeps(max_filter_rows)                                       &
                                   ! sweeps for 1-2-1 filter
     &, v_sweeps(max_filter_rows)                                       &
                                   ! sweeps for 1-2-1 filter
     &, horizontal_level                                                &
                               ! steep slope control
     &, diff_order_thermo                                               &
                           ! diffusion order for theta/w
     &, diff_order_wind                                                 &
                           ! diffusion order for u,v
     &, tar_horizontal                                                  &
                            ! steep slope control
     &, test_level                                                      &
                            ! test start level for targeted diffusion
     &, start_level                                                     &
                            ! start level for targeted diffusion
     &, end_level           ! end level for targeted diffusion

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_at_u (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:rows+halo_j, model_levels)                     &
     &, r_at_v (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:n_rows+halo_j, model_levels)

      Real :: exner_theta_levels                                        &
     &    (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)

      Integer                                                           &
     &  n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, proc_row_group                                                  &
     &, proc_col_group                                                  &
     &, neighbour(4)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN/OUT.
      Real                                                              &
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &       model_levels)                                              &
     &, R_w(row_length, rows, model_levels)                             &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
     &, theta_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &             model_levels)

!   Array arguments with intent(out):
      REAL, INTENT(OUT) ::                                              &
     & STASHwork(*)   ! Output array holding diagnostic fields

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

! local  variables
      Integer                                                           &
     &  i, j, k     ! Loop indices

! Local arrays
      Real                                                              &
     &  field( 1-off_x:row_length+off_x, 1-off_y:rows+off_y             &
     &                   , model_levels)

      Real :: T_inc (row_length, rows, model_levels)
      Real :: u_inc (row_length, rows, model_levels)
      Real :: v_inc (row_length, n_rows, model_levels)
      Real :: w_inc (row_length, rows, model_levels)
      Integer :: sect
      Integer :: item
      Integer :: im_index      !  internal model index for STASH arrays
      INTEGER :: Errorstatus = 0  ! initial value for error code
      CHARACTER (LEN=80) :: CMessage !  Error message

! ----------------------------------------------------------------------
! 0.1  Section for setting up diagnostic arrays
! ----------------------------------------------------------------------
      sect=13
      im_index    = internal_model_index(atmos_im)
      Cmessage    = ''

      IF (sf(0,sect)) THEN
        IF(sf(481,sect))t_inc(:,:,:)=theta_star(1:row_length, 1:rows, :)
        IF(sf(485,sect))u_inc(:,:,:)=R_u(1:row_length, 1:rows, :)
        IF(sf(486,sect))v_inc(:,:,:)=R_v(1:row_length,1:n_rows,:)
        IF(sf(487,sect))w_inc(:,:,:)=R_w(1:row_length, 1:rows, :)
      ENDIF

! ----------------------------------------------------------------------
! 1.0  Section for original polar filter
! ----------------------------------------------------------------------

      If (L_polar_filter_incs) Then

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('Polar_Filter_incs',3)

! DEPENDS ON: polar_filter_incs
        Call Polar_filter_incs(                                         &
     &                      theta, theta_star, R_u, R_v, R_w,           &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      north_lat_limit, south_lat_limit,           &
     &                      polar_filter_coefficient,                   &
     &                      polar_filter_n_sweeps, step_per_sweep,      &
     &                      polar_start_lat_limit,                      &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      n_proc, n_procx, n_procy, l_datastart,      &
     &                      neighbour, at_extremity, proc_row_group)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('Polar_Filter_incs',4)

      End If   ! L_polar_filter_incs

! ----------------------------------------------------------------------
! 2.0  Section for new polar filter
! ----------------------------------------------------------------------

      if( L_pftheta .or. L_diff_incs )then

          Do k = 1, model_levels
            Do j = 1-off_y, rows+off_y
              Do i = 1-off_x, row_length+off_x
                field(i,j,k) = theta_star(i,j,k) - theta(i,j,k)
              End Do
            End Do
          End Do !  k = 1, model_levels


! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   field, row_length, rows, model_levels,         &
     &                   off_x, off_y, fld_type_p, .false.)

        if( L_pofil_new )then

! DEPENDS ON: pfwthinc
          Call pfwthinc(                                                &
     &                  field, r_theta_levels, r_rho_levels,            &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  sec_theta_latitude, cos_v_latitude,             &
     &                  me, n_procy, delta_lambda, delta_phi,           &
     &                  rows, n_rows, row_length, model_levels,         &
     &                  max_filter_rows, u_begin, u_end,                &
     &                  u_sweeps, global_u_filter,                      &
     &                  model_levels, horizontal_level,                 &
     &                  diff_coeff_phi, diff_coeff_u, L_diff_incs )

        else  ! original combi filter

! DEPENDS ON: pofil_wth_incs
          Call pofil_wth_incs(                                          &
     &                      field, r_theta_levels, r_rho_levels,        &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      sec_theta_latitude, cos_v_latitude,         &
     &                      me, n_procy, delta_lambda, delta_phi,       &
     &                      rows, n_rows, row_length, model_levels,     &
     &                      max_filter_rows, u_begin, u_end,            &
     &                      u_sweeps, global_u_filter,                  &
     &                      model_levels, horizontal_level,             &
     &                      diff_order_thermo, diff_coeff_thermo,       &
     &                      diff_coeff_u, L_diff_incs )

        endif ! L_pofil_new

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                theta_star(i,j,k) = theta(i,j,k) + field(i,j,k)
              End Do
            End Do
          End Do

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     theta_star, row_length, rows, model_levels,  &
     &                     off_x, off_y, fld_type_p, .false.)

      endif ! L_pftheta .or. L_diff_incs

        if(L_pfw)then

          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                field(i,j,k) = R_w(i,j,k)
              End Do
            End Do
          End Do  ! k = 1, model_levels - 1

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   field, row_length, rows, model_levels - 1,     &
     &                   off_x, off_y, fld_type_p, .false.)

        if( L_pofil_new )then

! DEPENDS ON: pfwthinc
          Call pfwthinc(                                                &
     &                  field, r_theta_levels, r_rho_levels,            &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  sec_theta_latitude, cos_v_latitude,             &
     &                  me, n_procy, delta_lambda, delta_phi,           &
     &                  rows, n_rows, row_length, model_levels,         &
     &                  max_filter_rows, u_begin, u_end,                &
     &                  u_sweeps, global_u_filter,                      &
     &                  model_levels - 1, horizontal_level,             &
     &                  diff_coeff_phi, diff_coeff_u, L_diff_incs )

        else  ! original combi filter

! DEPENDS ON: pofil_wth_incs
          Call pofil_wth_incs(                                          &
     &                      field, r_theta_levels, r_rho_levels,        &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      sec_theta_latitude, cos_v_latitude,         &
     &                      me, n_procy, delta_lambda, delta_phi,       &
     &                      rows, n_rows, row_length, model_levels,     &
     &                      max_filter_rows, u_begin, u_end,            &
     &                      u_sweeps, global_u_filter,                  &
     &                      model_levels-1, horizontal_level,           &
     &                      diff_order_thermo, diff_coeff_thermo,       &
     &                      diff_coeff_u, L_diff_incs )

        endif ! L_pofil_new

          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                R_w(i,j,k) = field(i,j,k)
              End Do
            End Do
          End Do  ! k = 1, model_levels - 1

        endif ! L_pfw

      if( L_pfuv .or. L_diff_incs) then

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   R_u, row_length, rows, model_levels,           &
     &                   off_x, off_y, fld_type_u, .true.)

        if( L_pofil_new )then

! DEPENDS ON: pfuinc
          Call pfuinc(                                                  &
     &                R_u, r_at_u,                                      &
     &                r_theta_levels, r_rho_levels,                     &
     &                off_x, off_y, halo_i, halo_j,                     &
     &                sec_theta_latitude, cos_v_latitude,               &
     &                me, n_procy, delta_lambda, delta_phi,             &
     &                rows, n_rows, row_length, model_levels,           &
     &                max_filter_rows, u_begin, u_end, u_sweeps,        &
     &                global_u_filter, horizontal_level,                &
     &                diff_coeff_phi, diff_coeff_u, L_diff_incs )

        else  ! original combi filter

! DEPENDS ON: pofil_u_incs
          Call pofil_u_incs(                                            &
     &                 R_u, r_at_u,                                     &
     &                 r_theta_levels, r_rho_levels,                    &
     &                 off_x, off_y, halo_i, halo_j,                    &
     &                 sec_theta_latitude, cos_v_latitude,              &
     &                 me, n_procy, delta_lambda, delta_phi,            &
     &                 rows, n_rows, row_length, model_levels,          &
     &                 max_filter_rows, u_begin, u_end, u_sweeps,       &
     &                 global_u_filter,                                 &
     &                 horizontal_level, diff_order_wind,               &
     &                 diff_coeff_wind, diff_coeff_u, L_diff_incs )

        endif ! L_pofil_new

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   R_v, row_length, n_rows, model_levels,         &
     &                   off_x, off_y, fld_type_v, .true.)

        if( L_pofil_new )then

! DEPENDS ON: pfvinc
          Call pfvinc(                                                  &
     &                 R_v, r_at_v, r_theta_levels, r_rho_levels,       &
     &                 off_x, off_y, halo_i, halo_j,                    &
     &                 sec_v_latitude, cos_theta_latitude,              &
     &                 me, n_procy, at_extremity,                       &
     &                 delta_lambda, delta_phi,                         &
     &                 rows, n_rows, row_length, model_levels,          &
     &                 max_filter_rows, v_begin, v_end, v_sweeps,       &
     &                 global_v_filter,                                 &
     &                 horizontal_level, model_domain,                  &
     &                 diff_coeff_phi, diff_coeff_v, L_diff_incs )

        else  ! original combi filter

! DEPENDS ON: pofil_v_incs
          Call pofil_v_incs(                                            &
     &                 R_v, r_at_v, r_theta_levels, r_rho_levels,       &
     &                 off_x, off_y, halo_i, halo_j,                    &
     &                 sec_v_latitude, cos_theta_latitude,              &
     &                 me, n_procy, delta_lambda, delta_phi,            &
     &                 rows, n_rows, row_length, model_levels,          &
     &                 max_filter_rows, v_begin, v_end, v_sweeps,       &
     &                 global_v_filter,                                 &
     &                 horizontal_level, diff_order_wind,               &
     &                 diff_coeff_wind, diff_coeff_v, L_diff_incs )

        endif ! L_pofil_new

      endif !  L_pfuv .or. L_diff_incs

!*******************  increments for STASH   **************************

      IF (sf(0,sect)) THEN

! T increment
      item = 481
        IF(sf(item,sect)) THEN

          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                T_inc(i,j,k) = (theta_star(i,j,k) - T_inc(i,j,k))       &
     &                                  *exner_theta_levels(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        T_inc,                                                    &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! u wind increment
        item=485
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                u_inc(i,j,k) = R_u(i,j,k) - u_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        u_inc,                                                    &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! v wind increment
        item=486
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels
            DO j=1,n_rows
              DO i=1,row_length
                v_inc(i,j,k) = R_v(i,j,k) - v_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        v_inc,                                                    &
     &        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! w wind increment
        item=487
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                w_inc(i,j,k) = R_w(i,j,k) - w_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        w_inc,                                                    &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

      ENDIF ! sf(0,sect)

      return

      END SUBROUTINE NI_filter_incs_Ctl

