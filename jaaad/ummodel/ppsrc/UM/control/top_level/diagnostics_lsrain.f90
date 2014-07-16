
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      Subroutine diagnostics_lsrain(                                    &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels                           &
     &,                      lspice_dim1,lspice_dim2,lspice_dim3        &
     &,                      timestep                                   &
     &,                      at_extremity                               &
     &,                      L_DUST                                     &
     &,                      p_layer_centres                            &
     &,                      T, q, qcl, qcf, cf, cfl, cff               &
     &,                      T_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n &
     &,                      ls_rain, ls_snow                           &
     &,                      ls_rain3d,ls_snow3d,rainfrac3d             &
     &,                      RNOUT_TRACER,LSCAV_DUST_ALL,LSCAV_TR       &
     &,                      lscav_nh3                                  &
     &,                      rnout_soot, lscav_soot                     &
     &,                      rnout_bmass, lscav_bmass                   &
     &,                      rnout_ocff, lscav_ocff                     &
     &,                      l_mixing_ratio                             &
     &,                      PSDEP,PSAUT,PSACW,PSACR                    &
     &,                      PSACI,PSMLT,PSMLTEVP                       &
     &,                      PRAUT,PRACW,PREVP                          &
     &,                      PGAUT,PGACW,PGACS,PGMLT                    &
     &,                      PIFRW,PIPRM,PIDEP,PIACW                    &
     &,                      PIACR,PIMLT,PIMLTEVP                       &
     &,                      PIFALL,PSFALL,PRFALL,PGFALL                &
     &,                      PLSET,PLEVPSET                             &
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork                                                        &
     &  )

Use ac_diagnostics_mod, Only :  &
    LSRR, LSSR, TINC_PPN

! Purpose:
!          Calculates diagnostics generated from large scale
!          precipitation (UM section 4).
!
! Method:
!  Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the large scale
! precipitation routines called previously. Each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing, except where indicated.
! NOTE: Although moisture field diagnostics are available from this
! section (and correspond to equivalent variables calculated at UM4.5
! and earlier), the most appropriate place to extract moisture
! diagnostics is after section 9, which is the end of the moist
! processes calculations during the timestep.
!
!  Diagnostics currently available: (in order calculated)
!
! STASH item (all section 4 )
! ------------------------------------------------------------
! 181 temperature increment across ls precip routines (model levels)
! 182 humidity    increment across ls precip routines (model levels)
! 183 qcl         increment across ls precip routines (model levels)
! 184 qcf         increment across ls precip routines (model levels)
! 192 cf          increment across ls precip routines (model levels)
! 193 cfl         increment across ls precip routines (model levels)
! 194 cff         increment across ls precip routines (model levels)
! 201 large scale rain amount (kg/m2 per timestep)    (surface)
! 202 large scale snow amount (kg/m2 per timestep)    (surface)
! 203 large scale rainfall rate (kg/m2/s)             (surface)
! 204 large scale snowfall rate (kg/m2/s)             (surface)
!   4 temperature           after ls precip           (model levels)
! 205 cloud water (qcl)     after ls precip           (model levels)
! 206 cloud ice (qcf)       after ls precip           (model levels)
! 207 relative humidity     (percent)                 (model levels)
!  10 specific humidity (q) after ls precip           (model levels)
! 220 Large scale rainout of soot (kg/m2/s)           (surface)
! 221 Large scale washout of soot (kg/m2/s)           (surface)
! 223 snowfall rate                                   (model levels)
! 224 supercooled liquid water content                (model levels)
! 225 supercooled rainfall rate                       (model levels)
! 227 rain fraction                                   (model levels)
! 228 Large scale rainout of fossil-fuel organic carbon (kg/m2/s) (srf)
! 229 Large scale washout of fossil-fuel organic carbon (kg/m2/s) (srf)
! 237 Large scale rainout of biomass smoke (kg/m2/s)  (surface)
! 238 Large scale washout of biomass smoke (kg/m2/s)  (surface)
! Microphysical process rate diagnostics all on wet_model_levels
! 240 Homogeneous nucleation rate (kg/kg/s)
! 241 Heterogeneous nucleation rate (kg/kg/s)
! 243 Ice deposition rate (kg/kg/s)
! 245 Snow deposition rate (kg/kg/s)
! 247 Ice collection rate of cloud liquid water (riming) (kg/kg/s)
! 248 Snow collection rate of cloud liquid water (riming) (kg/kg/s)
! 249 Ice collection rate of rain (capture) (kg/kg/s)
! 250 Snow collection rate of rain (capture) (kg/kg/s)
! 251 Evaporation rate of melting ice (kg/kg/s)
! 252 Evaporation rate of melting snow (kg/kg/s)
! 253 Melting rate for ice (kg/kg/s)
! 254 Melting rate for snow (kg/kg/s)
! 255 Snow aggregate autoconversion rate (kg/kg/s)
! 256 Snow collection rate of ice (capture) (kg/kg/s)
! 257 Rain autoconversion rate (kg/kg/s)
! 258 Rain collection rate of cloud liquid water (accretion) (kg/kg/s)
! 259 Evaporation rate of rain (kg/kg/s)
! 260 Graupel autoconversion rate  (kg/kg/s)
! 261 Graupel collection rate of cloud water (accretion) (kg/kg/s)
! 262 Graupel collection rate of snow (capture) (kg/kg/s)
! 263 Melting rate for graupel (kg/kg/s)
! 265 Ice crystal sedimentation rate (kg/kg/s)
! 266 Snow sedimentation rate (kg/kg/s)
! 267 Rain sedimentation rate (kg/kg/s)
! 268 Graupel sedimentation rate (kg/kg/s)
! 269 Droplet settling rate of liquid (kg/kg/s)
! 270 Evaporation rate for settling droplets (kg/kg/s)
!
! History:
! Version   Date     Comment
! ----     -------     -------
! 5.0  30/11/99 Original version. J-C Thil.
! 5.1  09/12/99 Add error trapping. Rick Rawlins
! 05/04/00  5.1    Provide explicit model increments as STASH output
!                  diagnostics. R Rawlins
! 5.2  21/11/00 Add 3D precipitation field diagnostics. Damian Wilson
! 5.3  04/10/01 Add diagnostics for wet scavenged Sulphur Cycle
!               tracers (SO2 and SO4_DIS)         M Woodage, A Jones
! 5.3  14/08/01 Add (4,207) RH on model levels. R Rawlins
! 5.4  22/07/02 Add PC2 diagnostics.  Damian Wilson
! 5.4  10/07/02 Add (4,220) and (4,221) Rainout and washout of soot.
!                                                           P. Davison
! 5.4  05/09/02  Add diagnostic (4,215) for wet scavenged NH3 (for
!                      Sulphur Cycle).                       M Woodage
! 5.5  05/02/03 Add (4,237) and (4,238) Rainout and washout of biomass
!               smoke aerosol.                          P Davison
!
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
! 6.2  15/02/06 Fix to dust deposition diagnostic.  Stephanie Woodward
!
! 6.2  05/01/06 Add RH wrt water diagnostic (4,208). Richard Forbes.
! 6.2  31/01/06 Add microphysics process rate diags. Richard Forbes
! 6.4  14/08/06 Use mixing ratio version of qsat for 4207. D. Wilson
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      Implicit None
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


! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_DUST                                                          &
                        !switch for mineral dust
     &, l_mixing_ratio   ! Use mixing ratio formulation

! Parameters
!
      INTEGER, INTENT(IN) ::                                            &
     &  row_length                                                      &
                            ! number of points on a row
     &, rows                                                            &
                            ! number of rows in a theta field
     &, model_levels                                                    &
                            ! number of model levels
     &, wet_model_levels                                                &
                            ! number of model levels where moisture
     &, lspice_dim1                                                     &
                            ! Dimensions for 3D diagnostic arrays.
     &, lspice_dim2                                                     &
                            !  These are set to 1 in order to save
     &, lspice_dim3         !  memory if the diagnostics are not used.

      REAL, INTENT(IN) ::                                               &
     &  timestep


! Primary Arrays used in all models
      REAL, INTENT(IN) ::                                               &
     &  p_layer_centres(row_length, rows, 0:model_levels)               &
     &, T(row_length, rows, model_levels)                               &
     &, q(row_length,rows, wet_model_levels)                            &
     &, qcl(row_length, rows, wet_model_levels)                         &
     &, qcf(row_length, rows, wet_model_levels)                         &
     &, cf(row_length, rows, wet_model_levels)                          &
     &, cfl(row_length, rows, wet_model_levels)                         &
     &, cff(row_length, rows, wet_model_levels)                         &
! Time level n values for increment diagnostics
     &, T_n  (row_length, rows,     model_levels)                       &
     &, q_n  (row_length, rows, wet_model_levels)                       &
     &, qcl_n(row_length, rows, wet_model_levels)                       &
     &, qcf_n(row_length, rows, wet_model_levels)                       &
     &, cf_n  (row_length, rows, wet_model_levels)                      &
     &, cfl_n  (row_length, rows, wet_model_levels)                     &
     &, cff_n  (row_length, rows, wet_model_levels)                     &
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, ls_snow3d(lspice_dim1, lspice_dim2, lspice_dim3)                &
     &, rainfrac3d(lspice_dim1, lspice_dim2, lspice_dim3)
! Microphysics Process Rate diagnostic arrays
      REAL, Intent(In) ::                                               &
     &  PSDEP(row_length,rows,wet_model_levels)                         &
     &, PSAUT(row_length,rows,wet_model_levels)                         &
     &, PSACW(row_length,rows,wet_model_levels)                         &
     &, PSACR(row_length,rows,wet_model_levels)                         &
     &, PSACI(row_length,rows,wet_model_levels)                         &
     &, PSMLT(row_length,rows,wet_model_levels)                         &
     &, PSMLTEVP(row_length,rows,wet_model_levels)
      REAL, Intent(InOut) ::                                            &
     &  PRAUT(row_length,rows,wet_model_levels)                         &
     &, PRACW(row_length,rows,wet_model_levels)                         &
     &, PREVP(row_length,rows,wet_model_levels)
      REAL, Intent(InOut) ::                                            &
     &  PGAUT(row_length,rows,wet_model_levels)                         &
     &, PGACW(row_length,rows,wet_model_levels)                         &
     &, PGACS(row_length,rows,wet_model_levels)                         &
     &, PGMLT(row_length,rows,wet_model_levels)
      REAL, Intent(InOut) ::                                            &
     &  PIFRW(row_length,rows,wet_model_levels)                         &
     &, PIPRM(row_length,rows,wet_model_levels)                         &
     &, PIDEP(row_length,rows,wet_model_levels)                         &
     &, PIACW(row_length,rows,wet_model_levels)                         &
     &, PIACR(row_length,rows,wet_model_levels)                         &
     &, PIMLT(row_length,rows,wet_model_levels)                         &
     &, PIMLTEVP(row_length,rows,wet_model_levels)
      REAL, Intent(InOut) ::                                            &
     &  PIFALL(row_length,rows,wet_model_levels)                        &
     &, PSFALL(row_length,rows,wet_model_levels)                        &
     &, PRFALL(row_length,rows,wet_model_levels)                        &
     &, PGFALL(row_length,rows,wet_model_levels)
      REAL, Intent(InOut) ::                                            &
     &  PLSET(row_length,rows,wet_model_levels)                         &
     &, PLEVPSET(row_length,rows,wet_model_levels)


! Used as input and workspace
      REAL, INTENT(INOUT) ::                                            &
     &  ls_rain3d(lspice_dim1, lspice_dim2, lspice_dim3)                &
     &, rnout_tracer(lspice_dim1, lspice_dim2)                          &
     &, LSCAV_DUST_ALL(ROW_LENGTH,ROWS,NDIV)                            &
                                             !scavenged mineral dust

     &, lscav_tr(lspice_dim1, lspice_dim2)                              &
     &, lscav_nh3(lspice_dim1, lspice_dim2)                             &
     &, rnout_soot(lspice_dim1, lspice_dim2)                            &
     &, lscav_soot(lspice_dim1, lspice_dim2)                            &
     &, rnout_bmass(lspice_dim1, lspice_dim2)                           &
     &, lscav_bmass(lspice_dim1, lspice_dim2)                           &
     &, rnout_ocff(lspice_dim1, lspice_dim2)                            &
     &, lscav_ocff(lspice_dim1, lspice_dim2)

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
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------

! Diagnostic variables
      REAL, INTENT(INOUT) ::                                            &
     & STASHwork(*)     ! STASH workspace for section 4 (LS precip)

! Local variables
      Integer                                                           &
     & i, j, k, ji                                                      &
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      Integer sect,item    ! STASH section, item no.s
      Parameter (sect = 4) !  for microphysics - large scale rain

      Real work_3d(row_length, rows,model_levels) ! work space

      Character*80  cmessage

      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_lsrain')

      Integer                                                           &
     &  im_index        ! internal model index

! External routines
      External                                                          &
     & copydiag, copydiag_3d                                            &
     &  ,Ereport

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

!L Copy diagnostic information to STASHwork for STASH processing

! increment diagnostics= modified - previous

      item = 181  ! temperature increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = T(i,j,k) - T_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

! And set dry level increments to zero explicitly
         Do k = wet_model_levels+1,model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = 0.0
             End do ! i
           End do   ! j
         End do     ! k

        If (.not.allocated(TINC_PPN)) Then
          Allocate ( TINC_PPN(row_length*rows,model_levels) )
        End If

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 181)"//cmessage
         End if

         Do k = 1,model_levels
           Do j = 1,rows
             Do i = 1,row_length
               ji = (j-1)*row_length+i
               TINC_PPN(ji,k) = work_3d(i,j,k) 
             End do ! i
           End do   ! j
         End do     ! k

      End if  !  sf(item,sect)

      item = 182  ! humidity increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = q(i,j,k) - q_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

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

      End if  !  sf(item,sect)

      item = 183  ! qcl increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = qcl(i,j,k) - qcl_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

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

      End if  !  sf(item,sect)

      item = 184  ! qcf increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = qcf(i,j,k) - qcf_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 184)"//cmessage
         End if

      End if  !  sf(item,sect)

      item = 192  ! cf increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = cf(i,j,k) - cf_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

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

      End if  !  sf(item,sect)

      item = 193  ! cfl increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = cfl(i,j,k) - cfl_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

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

      End if  !  sf(item,sect)

      item = 194  ! cff increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = cff(i,j,k) - cff_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 194)"//cmessage
         End if

      End if  !  sf(item,sect)

! Item 201 Large scale rain

      If (icode <= 0 .and. sf(201,4)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(201,4,im_index)),ls_rain,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,201,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(ls_rain)"
         Endif
! Code to convert rate to amount for a given timestep

         Do i=1,row_length*rows
            STASHwork(si(201,4,im_index)+i-1)=                          &
     &           STASHwork(si(201,4,im_index)+i-1)* timestep
         End do

      End if


! Item 202 Large scale snow

      If (icode <= 0 .and. sf(202,4)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(202,4,im_index)),ls_snow,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,202,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(ls_snow)"
         End if

         Do i=1,row_length*rows
            STASHwork(si(202,4,im_index)+i-1)=                          &
     &           STASHwork(si(202,4,im_index)+i-1)* timestep
        End do

      End if


! Item 203 Large scale rain

      If (icode <= 0 .and. sf(203,4)) Then

      If (.not.allocated(LSRR)) Then
        Allocate ( LSRR(row_length*rows) )
      End If

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(203,4,im_index)),ls_rain,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,203,                                           &
     &       icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(ls_rain)"
         End if

         Do j = 1,rows
           Do i = 1,row_length
             ji = (j-1)*row_length+i
             LSRR(ji) = ls_rain(i,j)
           End Do
         End Do

      End if


! Item 204 Large scale snow

      If (icode <= 0 .and. sf(204,4)) then

      If (.not.allocated(LSSR)) Then
        Allocate ( LSSR(row_length*rows) )
      End If

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(204,4,im_index)),ls_snow,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,204,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(ls_snow)"
         End if

         Do j = 1,rows
           Do i = 1,row_length
             ji = (j-1)*row_length+i
             LSSR(ji) = ls_snow(i,j)
           End Do
         End Do

      End if



! Items 231-236 mineral dust scavenged by LS precip

      IF (L_DUST) THEN

        DO K = 1,NDIV

          IF (ICODE <= 0 .AND. SF(230+K,4)) THEN

!       Convert "per timestep" diagnostic to "per second":
            DO J=1,LSPICE_DIM2
              DO I=1,LSPICE_DIM1
                LSCAV_DUST_ALL(I, J, K)=LSCAV_DUST_ALL(I, J, K)/TIMESTEP
              ENDDO
            ENDDO

! DEPENDS ON: copydiag
            CALL COPYDIAG(STASHWORK(SI(230+K,4,IM_INDEX)),              &
     &         LSCAV_DUST_ALL(1:ROW_LENGTH,1:ROWS,K),                   &
     &         ROW_LENGTH,ROWS,0,0,0,0, AT_EXTREMITY,                   &
     &         ATMOS_IM,4,230+K,                                        &
     &         ICODE,CMESSAGE)

            IF (ICODE  >   0) THEN
              CMESSAGE=": ERROR IN COPYDIAG(LSCAV_DUST_ALL)"
            ENDIF

          ENDIF

        ENDDO !NDIV

      ENDIF !L_DUST


! Item 215 NH3 scavenged by LS precip

      If (icode <= 0 .and. sf(215,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_nh3(i, j)=lscav_nh3(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(215,4,im_index)),lscav_nh3,          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,215,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage="microphysics_ctl  : error in copydiag(lscav_nh3)"
        End if

      End if
!
! Item 216 SO2 scavenged by LS precip

      If (icode <= 0 .and. sf(216,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_tr(i, j)=lscav_tr(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(216,4,im_index)),lscav_tr,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,216,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(lscav_tr)"
        End if

      End if

! Item 219 Dissolved SO4 aerosol scavenged by LS precip

      If (icode <= 0 .and. sf(219,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            rnout_tracer(i, j)=rnout_tracer(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(219,4,im_index)),rnout_tracer,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,219,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(rnout_tracer)"
        End if

      End if
!
! Item 220 Soot scavenged by LS rainout

      If (icode <= 0 .and. sf(220,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            rnout_soot(i, j)=rnout_soot(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(220,4,im_index)),rnout_soot,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,220,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(rnout_soot)"
        End if

      End if

! Item 221 Soot scavenged by LS washout

      If (icode <= 0 .and. sf(221,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_soot(i, j)=lscav_soot(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(221,4,im_index)),lscav_soot,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,221,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(lscav_soot)"
        End if

      End if
!
! Item 228 Fossil-fuel organic carbon scavenged by LS rainout

      If (icode <= 0 .and. sf(228,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            rnout_ocff(i, j)=rnout_ocff(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(228,4,im_index)),rnout_ocff,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,228,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(rnout_ocff)"
        End if

      End if

! Item 229 Fossil-fuel organic carbon scavenged by LS washout

      If (icode <= 0 .and. sf(229,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_ocff(i, j)=lscav_ocff(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(229,4,im_index)),lscav_ocff,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,229,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(lscav_ocff)"
        End if

      End if
!
!L Copy T to STASHwork

      If (icode <= 0 .and. sf(004,4)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(004,4,im_index)),T,              &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,004,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,004,                                           &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(T)"
         End if

      End if



!L Copy Cloud water to STASHwork

      If (icode <= 0 .and. sf(205,4)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(205,4,im_index)),qcl,            &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,205,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,205,                                           &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(cloud water)"
         End if

      End if

      If (icode <= 0 .and. sf(206,4)) Then

!L Copy Cloud ice to STASHwork

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(206,4,im_index)),qcf,            &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,206,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,206,                                           &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(cloud ice)"
         End if

      End if

! -----------------------------
! 207 relative humidity wrt ice (T<0degC) and water (T>0degC) (mdl levs)

      item = 207
      If (icode <= 0 .and. sf(item,sect)) Then
! DEPENDS ON: qsat_mix
         Call qsat_mix(work_3d,T,p_layer_centres(1,1,1),                &
     &             row_length*rows*wet_model_levels,l_mixing_ratio)
         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = q(i,j,k)/work_3d(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
               If (work_3d(i,j,k)  <   0.0) Then
                 work_3d(i,j,k) = 0.0
               End if

             End do  ! i
           End do    ! j
         End do      ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 207)"//cmessage
         End if

      End if   ! item 207

! -------------------------------
! 208 relative humidity wrt water (on model levels)

      item = 208
      If (icode  <=  0 .and. sf(item,sect)) Then
         ! q saturation is put in work_3d array
! DEPENDS ON: qsat_wat
         Call qsat_wat(work_3d,T,p_layer_centres(1,1,1),                &
     &             row_length*rows*wet_model_levels )
         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = q(i,j,k)/work_3d(i,j,k)*100.
               !  Supersaturation wrt water is limited to =< 100%
               If (work_3d(i,j,k) > 100.0) Then
                 work_3d(i,j,k) = 100.0
               End if
               !  Negative humidity also removed from the diagnostic
               If (work_3d(i,j,k) < 0.0) Then
                 work_3d(i,j,k) = 0.0
               End if

             End do  ! i
           End do    ! j
         End do      ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d (item 208, rhw)"//cmessage
         End if

      End if   ! item 208
!L Copy q  to STASHwork

      If (icode <= 0 .and. sf(010,4)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(010,4,im_index)),q,              &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,010,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,010,                                           &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(q)"
         End if

      End if


! Copy ls_rain3d  to STASHwork
!
      If (icode <= 0 .and. sf(222,4)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(222,4,im_index)),ls_rain3d,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,222,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,222,                                           &
     &        icode,cmessage)
!
         If (icode >  0) Then
            cmessage=":error in copydiag_3d(ls_rain3d)"
         End if
!
      End if
!
!
! Copy ls_snoww3d  to STASHwork
!
      If (icode <= 0 .and. sf(223,4)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223,4,im_index)),ls_snow3d,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,223,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,223,                                           &
     &        icode,cmessage)
!
         If (icode >  0) Then
            cmessage=":error in copydiag_3d(ls_snow3d)"
         End if
!
      End if
!
!
! Need to produce diagnostic 225 before 224 in order to save memory.
!
      If(icode <= 0 .and. sf(225,4)) Then
!
!  Supercooled 3D rain content. It is equal to
!  the 3D rainrate at T < 0 and equal to 0 at T > 0
!  Alter the array LS_RAIN3D directly
!
          Do k=1,wet_model_levels
            Do j=1,rows
              Do i=1,row_length
                If (T(i,j,k)  >=  zerodegc) THEN
! Warm temperatures
                  ls_rain3d(i,j,k)=0.0
                End if
              End do
            End do
          End do
!
! Copy supercooled rain (now in ls_rain3d)  to STASHwork
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(225,4,im_index)),ls_rain3d,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,225,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,225,                                           &
     &        icode,cmessage)
!
         If (icode >  0) Then
            cmessage=":error in copydiag_3d(supercooled rain)"
         End if
!
      End if
!
      IF (icode <= 0 .and. SF(224,4)) THEN
!
!  Supercooled liquid water content at TIME LEVEL N. It is equal to
!  the liquid water content at T < 0 and equal to 0 at T > 0
!  Use LS_RAIN3D as the array in order to save memory
!
        Do k=1,wet_model_levels
          Do j=1,rows
            Do i=1,row_length
              If (T(i,j,k)  <   zerodegc) Then
! Supercooled temperatures
! Use time level n fields in this diagnostic
                ls_rain3d(i,j,k)=qcl_n(i,j,k)
              Else
! Warm temperatures
                ls_rain3d(i,j,k)=0.0
              End if
            End do
          End do
        End do
!
! Copy supercooled liquid (now in ls_rain3d)  to STASHwork
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(224,4,im_index)),ls_rain3d,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,224,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,224,                                           &
     &        icode,cmessage)
!
        If (icode >  0) Then
           cmessage=":error in copydiag_3d(supercooled liq)"
        End if
!
      End if
!
!
      If (icode <= 0 .and. sf(227,4)) Then
!
! Copy rain fraction to stashwork
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(227,4,im_index)),rainfrac3d,     &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,227,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,227,                                           &
     &        icode,cmessage)
!
         If (icode >  0) Then
            cmessage=":error in copydiag_3d(rain fraction)"
         End if
!
      End if
!
! Item 237 Biomass scavenged by LS rainout

      If (icode <= 0 .and. sf(237,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            rnout_bmass(i, j)=rnout_bmass(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(237,4,im_index)),rnout_bmass,        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,237,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(rnout_bmass)"
        End if

      End if

! Item 238 Biomass scavenged by LS washout

      If (icode <= 0 .and. sf(238,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_bmass(i, j)=lscav_bmass(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(238,4,im_index)),lscav_bmass,        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,238,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(lscav_bmass)"
        End if

      End if
!

        !---------------------------------------------------------------
        ! Homogeneous nucleation
        !---------------------------------------------------------------
        item = 240

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIFRW,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Heterogeneous nucleation
        !---------------------------------------------------------------
        item = 241

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIPRM,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Deposition of ice
        !---------------------------------------------------------------
        item = 243

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIDEP,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Deposition of snow aggregates
        !---------------------------------------------------------------
        item = 245

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSDEP,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If


        !---------------------------------------------------------------
        ! Ice collection of cloud liquid water (riming)
        !---------------------------------------------------------------
        item = 247

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIACW,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Snow collection of cloud liquid water (riming)
        !---------------------------------------------------------------
        item = 248

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSACW,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Ice collection of rain (capture)
        !---------------------------------------------------------------
        item = 249

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIACR,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Snow collection of rain (capture)
        !---------------------------------------------------------------
        item = 250

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSACR,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If
        !---------------------------------------------------------------
        ! Evaporation of melting ice
        !---------------------------------------------------------------
        item = 251

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIMLTEVP,                                                 &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Evaporation of melting snow
        !---------------------------------------------------------------
        item = 252

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSMLTEVP,                                                 &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Melting ice
        !---------------------------------------------------------------
        item = 253

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIMLT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Melting snow
        !---------------------------------------------------------------
        item = 254

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSMLT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Snow aggregate autoconversion
        !---------------------------------------------------------------
        item = 255

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSAUT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Snow collection of ice (capture)
        !---------------------------------------------------------------
        item = 256

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSACI,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Rain autoconversion
        !---------------------------------------------------------------
        item = 257

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PRAUT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Rain collection of cloud liquid water (accretion)
        !---------------------------------------------------------------
        item = 258

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PRACW,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Evaporation of rain
        !---------------------------------------------------------------
        item = 259

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PREVP,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Graupel autoconversion
        !---------------------------------------------------------------
        item = 260

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PGAUT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Graupel collection of cloud liquid water (accretion)
        !---------------------------------------------------------------
        item = 261

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PGACW,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Graupel collection of snow (capture)
        !---------------------------------------------------------------
        item =262

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PGACS,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Melting graupel
        !---------------------------------------------------------------
        item = 263

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PGMLT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Ice crystal sedimentation
        !---------------------------------------------------------------
        item = 265

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIFALL,                                                   &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Snow sedimentation
        !---------------------------------------------------------------
        item = 266

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSFALL,                                                   &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Rain sedimentation
        !---------------------------------------------------------------
        item = 267

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PRFALL,                                                   &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Graupel sedimentation
        !---------------------------------------------------------------
        item = 268

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PGFALL,                                                   &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Droplet settling of liquid
        !---------------------------------------------------------------
        item = 269

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PLSET,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Evaporated settled droplets
        !---------------------------------------------------------------
        item = 270

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PLEVPSET,                                                 &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

! Single point exception handling
      If (icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      End if

      Return
      END SUBROUTINE diagnostics_lsrain

