
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine Interface:
      Subroutine Oper_Emergency

      IMPLICIT NONE
!
!    Method:  Reset UMUI determined namelist items to alternate
!             predetermined values dependent on value of Environment
!             variables known to the operational suite
!
! Current Code Owner:  Tim Westmacott
!
! History:
! version  date      comment
! -------  ----      -------
!  4.5    3/09/98   New.  Stuart Bell
!  5.0   21/06/99   Amendments for 'C-P C dynamics' grid upgrade; wind
!                   limit test no longer supported. R Rawlins
!  5.1   13/04/00   New (temporary) RUN_ASSIM_MODE for FASTRUNs.
!                   Adam Clayton.
!  5.1   22/02/00   Move PARVARS for TYPSIZE                 P.Burton
!  5.3   29/11/00   Define run length changes wrt analysis time for
!  5.3              global runs. Adam Clayton
!  5.4   12/08/02   ONLYTO48 added. ONLYTO6 renamed ONLYTO9. A Clayton
!  5.4   22/03/02   Remove comment on same line as #include
!                                                S. Carroll
!  6.2   05/01/06   Move ONLYTO9 outside #if defined(global)
!                   pre-processor directive M.Bush
!  6.2   05/01/06   Enable ONLYTO54 switch for Mesoscale ALTERNATE runs
!                   J.Bornemann

!
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: Control
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
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
! Print status information
! CPRINTST defines print status for standard output messages

      ! Minimum output, only essential messages
      INTEGER,PARAMETER :: PrStatus_Min    = 1

      ! Normal informative messages + warnings
      INTEGER,PARAMETER :: PrStatus_Normal = 2

      ! Operational status, all informative messages
      INTEGER,PARAMETER :: PrStatus_Oper   = 3

      ! All informative + extra diagnostic messages
      INTEGER,PARAMETER :: PrStatus_Diag   = 4

      INTEGER PrintStatus ! Control volume of standard output messages
      COMMON/PrintSt/PrintStatus

! CPRINTST end

! Local scalars:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='OperatorEmergency')
      INTEGER ICODE            ! return code from FORT_GET_ENV
      CHARACTER*80  ONLYTO3    ! value of EV MES SHORT RUN
      CHARACTER*80  ONLYTO9    ! value of EV SHORT RUN
      CHARACTER*80  ONLYTO12   ! value of EV MOGREPS SHORT RUN
      CHARACTER*80  ONLYTO15   ! value of EV MOGREPS SHORT RUN
      CHARACTER*80  ONLYTO18   ! value of EV MOGREPS SHORT RUN
      CHARACTER*80  ONLYTO21   ! value of EV MOGREPS SHORT RUN
      CHARACTER*80  ONLYTO24   ! value of EV FOAM SHORT RUN
      CHARACTER*80  ONLYTO42   ! value of EV EXTENDED UK4
      CHARACTER*80  ONLYTO48   ! value of EV GL 2DAY RUN
      CHARACTER*80  ONLYTO54   ! value of EV MES ALTERNATE RUN
      CHARACTER*80  ONLYTO72   ! value of EV GL 3DAY RUN
      CHARACTER*80  FASTRUN    ! value of EV FASTRUN
      CHARACTER*80  SHORTSTEP  ! value of EV SHORTSTEP

! Function & Subroutine calls:
      External FORT_GET_ENV

! FASTRUN
       CALL FORT_GET_ENV('FASTRUN',7,FASTRUN,80,ICODE)
       if(FASTRUN == 'true'.AND.ICODE == 0)then
        RUN_ASSIM_MODE='FastRun   '
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &    "FASTRUN=true, setting RUN_ASSIM_MODE=",RUN_ASSIM_MODE
       ENDIF  ! PrintStatus test
       endif

! SHORTSTEP
       CALL FORT_GET_ENV('SHORTSTEP',9,SHORTSTEP,80,ICODE)
       if(SHORTSTEP == 'true'.AND.ICODE == 0)then
        IF(PrintStatus >= PrStatus_Oper) THEN
         WRITE(6,*) RoutineName,':Warning, Wind limit re-setting '      &
     &             ,'for SHORTSTEP no longer supported'
        ENDIF  ! PrintStatus test
       endif

! UM6.5 : MODEL_ANALYSIS_HRS replaced by MODEL_ANALYSIS_MINS
! ONLYTO72 for GL atmos plus 6 hour assm
       CALL FORT_GET_ENV('ONLYTO72',8,ONLYTO72,80,ICODE)
       if(ONLYTO72 == 'true'.AND.ICODE == 0)then
       RUN_TARGET_END(3)=3                    ! days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0 +3 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*72)-           &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &     "ONLYTO72=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif

! ONLYTO48 for GL atmos
       CALL FORT_GET_ENV('ONLYTO48',8,ONLYTO48,80,ICODE)
       if(ONLYTO48 == 'true'.AND.ICODE == 0)then
       RUN_TARGET_END(3)=2                  ! days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*48)-           &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &     "ONLYTO48=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif

! ONLYTO21 for MOGREPS GL atmos
       CALL FORT_GET_ENV('ONLYTO21',8,ONLYTO21,80,ICODE)
       if(ONLYTO21 == 'true'.AND.ICODE == 0)then
       RUN_TARGET_END(3)=0                  ! days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+21 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*21)-           &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &     "ONLYTO21=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif

! ONLYTO18 for MOGREPS GL atmos
       CALL FORT_GET_ENV('ONLYTO18',8,ONLYTO18,80,ICODE)
       if(ONLYTO18 == 'true'.AND.ICODE == 0)then
       RUN_TARGET_END(3)=0                  ! days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+18 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*18)-           &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &     "ONLYTO18=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif

!  ONLYTO9 for atmos
      CALL FORT_GET_ENV('ONLYTO9',7,ONLYTO9,80,ICODE)
      if(ONLYTO9 == 'true'.AND.ICODE == 0)then
      RUN_TARGET_END(3)=0                    ! days
      RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+9 ! hours
      RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*9)-             &
     &                  60.0*RUN_TARGET_END(4)           ! minutes 
      IF(PrintStatus >= PrStatus_Oper) THEN
      if(mype == 0)WRITE(6,*)                                           &
     &   "ONLYTO9=true, setting RUN_TARGET_END=",RUN_TARGET_END
      ENDIF  ! PrintStatus test
      endif
!  ONLYTO12 for atmos
      CALL FORT_GET_ENV('ONLYTO12',8,ONLYTO12,80,ICODE)
      if(ONLYTO12 == 'true'.AND.ICODE == 0)then
      RUN_TARGET_END(3)=0                    ! days
      RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+12 ! hours
      RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*12)-            &
     &                  60.0*RUN_TARGET_END(4)            ! minutes 
      IF(PrintStatus >= PrStatus_Oper) THEN
      if(mype == 0)WRITE(6,*)                                           &
     &   "ONLYTO12=true, setting RUN_TARGET_END=",RUN_TARGET_END
      ENDIF  ! PrintStatus test
      endif
!  ONLYTO15 for atmos
      CALL FORT_GET_ENV('ONLYTO15',8,ONLYTO15,80,ICODE)
      if(ONLYTO15 == 'true'.AND.ICODE == 0)then
      RUN_TARGET_END(3)=0                    ! days
      RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+15 ! hours
      RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*15)-            &
     &                  60.0*RUN_TARGET_END(4)            ! minutes 
      IF(PrintStatus >= PrStatus_Oper) THEN
      if(mype == 0)WRITE(6,*)                                           &
     &   "ONLYTO15=true, setting RUN_TARGET_END=",RUN_TARGET_END
      ENDIF  ! PrintStatus test
      endif

      RETURN
      END SUBROUTINE Oper_Emergency