
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Routine : TEMPHIST -------------------------------------------------
!LL
!LL Purpose :Write current contents of history common block to temporary
!LL          or interim history file - overwriting previous record
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  Author:   A.B.SANGSTER       Date:           20 January 1990
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL   3.1  05/02/93    Portable Fortran unit no assigns
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL  3.5  06/04/95  Sub-Models stage 1: revise History and Control file
!LL                 contents.  RTHBarnes.
!LL  4.5  05/05/98  Add DELIM='APOSTROPHE' to OPEN statement so that
!LL                 History file is written in correct form from
!LL                 Fujitsu.                       RBarnes@ecmwf.int
!    5.4  09/04/02   Need to specify DELIM='APOSTROPHE' in OPEN for
!                    namelists for Linux compilers
!                                                        S. Carroll
!    5.5  20/02/03  Modified to include NEC cpp        A.A.Dickinson
!    6.0  11/09/03  Modified to include IBM and SGI pre-processor
!                   directives (supplied by Zoe Chaplin)   P.Dando
!    6.1  18/08/04  Add DELIM=APOSTROPHE in open statement for ALTIX
!                   (supplied by Z. Chaplin).              P.Dando
!    6.2  16/05/05  Add DELIM=APOSTROPHE in open statement for HP.
!                   Simon Wilson
!    6.2  31/01/06   Changed defs from specific compilers to a general
!                    GNU/Linux option. T. Edwards
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered: H3,H40
!LL
!LL  Project task: H
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LL
!*L  Interface and arguments:
!
      SUBROUTINE TEMPHIST                                               &
     &         ( UNITHIST,ICODE,CMESSAGE )
!
      IMPLICIT NONE
!
      INTEGER       UNITHIST ! In  - Temporary history file unit
      INTEGER       ICODE    ! Out - Return code from routine
      CHARACTER*80  CMESSAGE ! Out - Return message if failure occured
!*
!
!L Common blocks
!
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
!
       CHARACTER *80 FILENAME
       CHARACTER *8  NLNAME
!*L EXTERNAL subroutines called
      EXTERNAL GET_FILE
!*
!
!L
!L 1. Open, rewind and write a new record
!L
      CALL GET_FILE(UNITHIST,FILENAME,80,ICODE)
       OPEN(UNITHIST,FILE=FILENAME,FORM='FORMATTED',                    &
     &  DELIM='APOSTROPHE',IOSTAT=ICODE)
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='TEMPHIST: Failed in OPEN of history file'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
        WRITE(6,*)'TEMPHIST: Warning message on OPEN of history file'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      REWIND(UNITHIST)

      NLNAME = 'NLIHISTO'
      WRITE(UNITHIST,NLIHISTO,ERR=200)
      NLNAME = 'NLCHISTO'
      WRITE(UNITHIST,NLCHISTO,ERR=200)
      NLNAME = 'NLIHISTG'
      WRITE(UNITHIST,NLIHISTG,ERR=200)
      NLNAME = 'NLCHISTG'
      WRITE(UNITHIST,NLCHISTG,ERR=200)
      NLNAME = 'NLCFILES'
      WRITE(UNITHIST,NLCFILES,ERR=200)

      go to 999
!
! Check for error
!
! Write error
  200 continue
      ICODE = 2
      CMESSAGE='TEMPHIST: Write ERROR on history file for namelist'//   &
     & NLNAME
!
 999  CONTINUE
!L
!L 2. Close and return
!L
      CLOSE(UNITHIST)
      RETURN
      END SUBROUTINE TEMPHIST
