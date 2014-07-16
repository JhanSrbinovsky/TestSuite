
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: GET_NAME -------------------------------------------------
!LL
!LL  Purpose: Generates an output file name of up to 14 characters using
!LL           the defined file naming convention, taking account of
!LL           file type, validity time, etc.
!LL           Obeys new filenaming convention introduced at version 2.7.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Author:   R A Stratton
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  27/01/93 : correct error in 2 character months - change jan
!LL                   from jn to ja so that no clash with june.
!LL  3.1 2/02/93 : added comdeck CHSUNITS to define NUNITS for i/o
!LL
!LL   3.1  22/02/93 : Cater for filename changes for boundary datasets.
!LL                   FILETYPE=u-z for LAM areas 1-6. D. Robinson
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.2  15/04/93  Correct Y_HUNDS in Absolute_long convention (TCJ).
!LL   3.2  08/07/93  Correct 1 T/S offset in reinit'ed pp names  (TCJ).
!LL   3.3  08/02/94  Modify calls to TIME2SEC/SEC2TIME to output/input
!LL                  elapsed times in days & secs, for portability. TCJ
!LL   3.4  17/06/94  Argument LCAL360 added and passed to SEC2TIME
!LL                                                        S.J.Swarbrick
!LL  4.1  30/07/96  Introduce Wave sub-model.  M Holt
!LL  4.4  11/07/97  Allow character filenames for PP files
!LL                 reinitialised on real month boundaries.  M Gallani
!LL  4.5  29/07/98  New naming convention for reinitialised boundary
!LL                 files. D. Robinson.
!LL  5.2  13/10/00  Dumpname fixed when MODEL_STATUS=SCS. Stuart Bell
!LL  5.2  16/01/01  New naming convention for reinitialized Macro
!LL                 files. (R.Hatcher)
!LL  5.2  19/07/99  Correct new filenaming convention for reinitialised
!LL                 boundary files. D. Robinson.
!LL  6.2  27/10/05  Add "Sub-hourly" filename convention. R Barnes.
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: S51
!LL
!LL  Project task: S51
!LL
!LL  External documentation: UM documentation paper 7 - Filenaming
!LL                          conventions for the Unified Model
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE GET_NAME(EXPT_ID,JOB_ID,ISUBMODEL,MEANLEV,TOGGLE,      &
     &          REINIT_STEPS,FILETYPE,LETTER_3, MODEL_STATUS,           &
     &     TIME_CONVENTION,ANALYSIS_HRS,FILENAME,ICODE,CMESSAGE,        &
     &     LCAL360)
!
      IMPLICIT NONE
      LOGICAL LCAL360
!
      CHARACTER*4   EXPT_ID     ! IN  - Experiment ident or alias
      CHARACTER*1   JOB_ID      ! IN  - Job ident within experiment
      INTEGER       ISUBMODEL   ! IN  - Submodel indicator
      INTEGER       MEANLEV     ! IN  - Mean level indicator
      INTEGER       TOGGLE      ! IN  - Alternately 1/2 for partial sums
      REAL          ANALYSIS_HRS! IN  - Hrs from basis time to analysis
                                ! UM6.5 - Allow for fractional ANALYSIS_HRS - 
                                !          change from INTEGER to REAL
      INTEGER       REINIT_STEPS! IN  - timesteps between file reinit-
!                                      ialisation for non-mean pp files
!                                      or -ve for Gregorian reinit.

      CHARACTER*1   FILETYPE    ! IN  - Code for file type
      CHARACTER*1   LETTER_3    ! IN  - character for use in position 9
!                                       of non-mean pp files.
      CHARACTER*14  MODEL_STATUS! IN  - Operational/NonOperational
!
      CHARACTER*17  TIME_CONVENTION ! IN  - Relative/Timestep/
!                                       Absolute_standard/Absolute_long/
!                                       Absolute_short
      CHARACTER*14  FILENAME    ! OUT - Generated file name
      INTEGER ICODE             ! OUT - Error return code
      CHARACTER*80 CMESSAGE
!
!*----------------------------------------------------------------------
!  Common blocks
!
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
!
! External subroutines called
!
      EXTERNAL SEC2TIME,DAY2CHAR,STP2TIME
!
!  Local variables
!
      INTEGER                                                           &
     & YYYY,MM,DD,HH,IMIN,ISEC   ! Current time values for filename
      INTEGER     COUNT,                                                &
                                 ! Counter for steps or hours
     &            DAYNO,                                                &
                                 ! day number
     &            DAYS,                                                 &
                                 ! Number of days for period
     &            HOURS,                                                &
                                 ! Number of hours for period
     &            I,                                                    &
                                 ! loop counter
     &            STEPS          ! number of steps
      INTEGER     END_DAYS       ! number of whole days from run start
      INTEGER     END_SECONDS    ! number of extra secs from run start
      INTEGER     MON            ! month for mean period
      INTEGER     A_STEPS_PER_HR    ! steps per hour for atmos sub-model
!
      CHARACTER*2 QW             ! Operational file prefix
!
      CHARACTER*1 SUBMODEL_ID    ! character for model a or o
      CHARACTER*1 FILETYPE_2     ! letter after FILETYPE in name
      CHARACTER*1 MEAN_PERIOD(4) ! default letter for mean period
!
      CHARACTER*1 Y_HUNDS        ! Character year identifier (hundreds)
      CHARACTER*1 Y_TENS         ! Character year identifier (tens)
      CHARACTER*1 Y_UNITS        ! Character year identifier (units)
      CHARACTER*1 M              ! Character month identifier
      CHARACTER*1 D              ! Character day-of-month identifier
      CHARACTER*1 H              ! Character hour identifier
      CHARACTER*1 HUNDREDS       ! Character for hundreds counter
      CHARACTER*1 TENS           ! Character for tens counter
      CHARACTER*1 UNITS          ! Character for units counter
      CHARACTER*1 DECI           ! Character for tenths (=mins)
      CHARACTER*1 CHAR_ID(36)    ! Valid characters for above (lookup)
      CHARACTER*1 SEPARATOR      ! character used as separator in name
      CHARACTER*1 STYLE          ! style of date in filename
      CHARACTER*3 CDAYNO         ! character day number
      CHARACTER*3 MONTH_3CHAR(12)! 3 character month identifier
      CHARACTER*2 MONTH_2CHAR(12)! 2 character month identifier
      CHARACTER*3 SEASON_3CHAR(12)! 3 character season identifier
      CHARACTER*2 SEASON_2CHAR(12)! 2 character season identifier
!
      DATA QW / 'qw'/
      DATA MEAN_PERIOD / '1', '2', '3', '4' /
      DATA CHAR_ID/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',   &
     &              'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',   &
     &              'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',   &
     &              'u', 'v', 'w', 'x', 'y', 'z' /
      DATA MONTH_3CHAR/ 'jan', 'feb', 'mar', 'apr', 'may', 'jun',       &
     &                  'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/
      DATA MONTH_2CHAR/ 'ja', 'fb', 'mr', 'ar', 'my', 'jn',             &
     &                  'jl', 'ag', 'sp', 'ot', 'nv', 'dc'/
      DATA SEASON_3CHAR/ 'ndj', 'djf', 'jfm', 'fma', 'mam', 'amj',      &
     &                   'mjj', 'jja', 'jas', 'aso', 'son', 'ond'/
      DATA SEASON_2CHAR/ 'nj', 'df', 'jm', 'fa', 'mm', 'aj',            &
     &                   'mj', 'ja', 'js', 'ao', 'sn', 'od'/
!L
!L----------------------------------------------------------------------
!L 1. Determine submodel id - (used in filechar 6 or 7 if operational)
!L
      IF (ISUBMODEL == ATMOS_SM) THEN
       SUBMODEL_ID= 'a'
      ELSE IF (ISUBMODEL == OCEAN_SM) THEN
       SUBMODEL_ID= 'o'
      ELSE IF (ISUBMODEL == WAVE_SM) THEN
       SUBMODEL_ID= 'w'
      ELSE
       ICODE=2
       CMESSAGE='GET_NAME: Illegal sub-model specified'
       GOTO 999
      ENDIF
      IF (ISUBMODEL == ATMOS_SM) THEN
! 1.1 Compute steps per hour for atmosphere sub_model
      A_STEPS_PER_HR = 3600*STEPS_PER_PERIODim(a_im)/                   &
     &                       SECS_PER_PERIODim(a_im)
      ENDIF
!
!L----------------------------------------------------------------------
!L 2. Determine style filename and separator
!L
      IF (FILETYPE /= 's') THEN
!L
!L 2.1 Relative time convention
!L
        IF (TIME_CONVENTION == 'Relative        ') THEN
          SEPARATOR='_'
          STYLE='A'
          IF (ISUBMODEL == atmos_sm) THEN ! atmosphere
          COUNT = STEPim(a_im) / A_STEPS_PER_HR - ANALYSIS_HRS
          ELSE IF (ISUBMODEL == ocean_sm) THEN ! ocean
            COUNT = STEPim(o_im) * SECS_PER_PERIODim(o_im) /            &
     &                 STEPS_PER_PERIODim(o_im) / 3600                  &
     &             -ANALYSIS_HRS
          ELSE IF (ISUBMODEL == wave_sm) THEN ! WAVE
            COUNT = STEPim(w_im) * SECS_PER_PERIODim(w_im) /            &
     &                 STEPS_PER_PERIODim(w_im) / 3600                  &
     &             -ANALYSIS_HRS
          ENDIF
!L
!L 2.2 Step time convention
!L
        ELSE IF (TIME_CONVENTION == 'Timestep         ') THEN

          SEPARATOR='_'
          STYLE='A'
          IF (ISUBMODEL == atmos_sm) THEN ! Atmosphere
            COUNT = STEPim(a_im)
          ELSE if (ISUBMODEL == ocean_sm) then ! Ocean
            COUNT = STEPim(o_im)
          ELSE if (ISUBMODEL == wave_sm) then ! WAVE
            COUNT = STEPim(w_im)
          ENDIF
!L
!L 2.3 Absolute time convention -standard version
!L
        ELSE IF (TIME_CONVENTION == 'Absolute_standard') THEN
          SEPARATOR='.'
          STYLE='B'
!L
!L 2.4 Absolute time convention - short
!L
        ELSE IF (TIME_CONVENTION == 'Absolute_short   ') THEN
          SEPARATOR='-'
          STYLE='B'
!L
!L 2.5 Absolute time convention - long
!L
        ELSE IF (TIME_CONVENTION == 'Absolute_long    ') THEN
          SEPARATOR='@'
          STYLE='B'
!L
!L 2.6 Sub-hourly filenaming time convention
!L
        ELSE IF (TIME_CONVENTION == 'Sub-hourly      ') THEN
          SEPARATOR='_'
          STYLE='A'
          IF (ISUBMODEL == atmos_sm) THEN ! atmosphere
            COUNT = STEPim(a_im) * SECS_PER_STEPim(a_im)
      write(6,*)'GET_NAM - COUNT ',a_im,COUNT,STEPim(a_im),             &
     & SECS_PER_STEPim(a_im)
          ELSE IF (ISUBMODEL == ocean_sm) THEN ! ocean
            COUNT = STEPim(o_im) * SECS_PER_STEPim(o_im)
          ELSE IF (ISUBMODEL == wave_sm) THEN ! WAVE
            COUNT = STEPim(w_im) * SECS_PER_STEPim(w_im)
          ENDIF
!
!   UM6.5 - ANALYSIS_HRS changed to REAL - 
!                  requires recoding of sub-hourly file-naming.  
             COUNT = 100*(REAL(COUNT)/3600.0 - ANALYSIS_HRS)
             COUNT = (COUNT/100)*100 + (REAL(MOD(COUNT,100))/100.0)*60

      write(6,*)'GET_NAM - hhmm ',COUNT,ANALYSIS_HRS
        ELSE
            ICODE=1
            CMESSAGE='GET_NAME: Illegal TIME_CONVENTION specified'
            GOTO 999
        ENDIF
!L----------------------------------------------------------------------
!L
!L 3.0 work out encoding of date time and filetype_2
!L
        IF (STYLE == 'A') THEN
          IF (COUNT <  0) THEN
            COUNT=-COUNT
            FILETYPE_2='z'
          IF (TIME_CONVENTION == 'Sub-hourly      ') THEN
!!!            FILETYPE_2='h'
            IF (COUNT  >=  36000) THEN
              ICODE = COUNT
              CMESSAGE='GET_NAME: COUNT too big for '//                 &
     &                 'sub-hourly filenaming convention'
              GOTO 999
            END IF
            HUNDREDS = CHAR_ID(MOD(COUNT/1000,36)+1)
            TENS =     CHAR_ID(MOD(COUNT/100, 10)+1)
            UNITS =    CHAR_ID(MOD(COUNT/10,  10)+1)
            DECI =     CHAR_ID(MOD(COUNT,     10)+1)
          ELSE
            IF (COUNT  >=  3600) THEN
              ICODE = COUNT
              CMESSAGE='GET_NAME: COUNT too big for '//                 &
     &                 'hourly or timestep filenaming convention'
              GOTO 999
            END IF
            HUNDREDS = CHAR_ID(MOD(COUNT/100,36)+1)
            TENS =     CHAR_ID(MOD(COUNT/10 ,10)+1)
            UNITS =    CHAR_ID(MOD(COUNT,    10)+1)
          END IF
          ELSE
            IF (FILETYPE == 'p') THEN
              FILETYPE_2=LETTER_3
        ELSE IF (FILETYPE == 'b') THEN   !  Boundary File
          FILETYPE_2=LETTER_3
            ELSE IF (FILETYPE == 'c') THEN   ! Macro File
              FILETYPE_2=LETTER_3
            ELSE
              FILETYPE_2='a'
            ENDIF
          IF (TIME_CONVENTION == 'Sub-hourly      ') THEN
!!!            FILETYPE_2='h'
            IF (COUNT  >=  36000) THEN
              ICODE = COUNT
              CMESSAGE='GET_NAME: COUNT too big for '//                 &
     &                 'sub-hourly filenaming convention'
              GOTO 999
            END IF
            HUNDREDS = CHAR_ID(MOD(COUNT/1000,36)+1)
            TENS =     CHAR_ID(MOD(COUNT/100, 10)+1)
            UNITS =    CHAR_ID(MOD(COUNT/10,  10)+1)
            DECI =     CHAR_ID(MOD(COUNT,     10)+1)
          ELSE
            IF (COUNT  >=  3600) THEN
              ICODE = COUNT
              CMESSAGE='GET_NAME: COUNT too big for '//                 &
     &                 'hourly or timestep filenaming convention'
              GOTO 999
            END IF
            HUNDREDS = CHAR_ID(MOD(COUNT/100,36)+1)
            TENS =     CHAR_ID(MOD(COUNT/10 ,10)+1)
            UNITS =    CHAR_ID(MOD(COUNT,    10)+1)
          END IF
          ENDIF
        ELSE IF (STYLE == 'B') THEN   ! some sort of absolute time
!
! Current date time is
!
          YYYY=I_YEAR
          MM  =I_MONTH
          DD  =I_DAY
          HH  =I_HOUR
          DAYNO = I_DAY_NUMBER
!
! Instantaneous files
!
          IF (MEANLEV == 0) THEN
!  dumps
            IF (FILETYPE == 'd') THEN
              FILETYPE_2 = 'a'
            ENDIF
!
! Work out reintialisation period for pp and boundary files.
! Note assumes reinitialisation period is whole number of hours. This
! is not strictly true but is probably ok for this purpose.
!
         IF (  FILETYPE == 'p'                                          &
                                 !  PP File
     &    .or. FILETYPE == 'b'                                          &
                                 !  Boundary File
     &    .or. FILETYPE == 'c'                                          &
                                 !  Macro File
     & ) THEN
              IF (ISUBMODEL == atmos_sm) then
                HOURS = REINIT_STEPS/A_STEPS_PER_HR
              ELSE IF (ISUBMODEL == ocean_sm) THEN
                HOURS = REINIT_STEPS*SECS_PER_PERIODim(o_im)            &
     &                             /(STEPS_PER_PERIODim(o_im)*3600)
              ELSE IF (ISUBMODEL == wave_sm) THEN
                HOURS = REINIT_STEPS*SECS_PER_PERIODim(w_im)            &
     &                             /(STEPS_PER_PERIODim(w_im)*3600)
              ENDIF
              if (REINIT_STEPS <  0) then ! Gregorian reinitialisation
                HOURS=720 ! dummy: could be anything divisible by 24
              endif
!   do further checks if multiple of 1 day
              IF (MOD(HOURS,24) == 0) THEN  ! whole days in period, or
                DAYS=HOURS/24               ! Gregorian reinit.
                IF (FILETYPE == 'b'                                     &
                                      !  Boundary File
     &           .or. FILETYPE == 'c'                                   &
                                      !  Macro File
     &        ) THEN
                  FILETYPE_2=LETTER_3
                ENDIF
                IF (DD == 1 .and. TIME_CONVENTION /= 'Absolute_short   '&
     &              .and. (DAYS == 30 .or. REINIT_STEPS <  0)) then
!  Original code didn't allow style=C for 3-month reinit files but new
!  code (Gregorian reinit) does, at least in section 3.0.
                  STYLE='C'                ! month in characters
                ENDIF
              ELSE
                IF (FILETYPE == 'b'                                     &
                                      !  Boundary File
     &           .or. FILETYPE == 'c'                                   &
                                      !  Macro File
     &        ) THEN
                  FILETYPE_2=LETTER_3
                ENDIF
              ENDIF
!
! For instantaneous pp file need to work out end of period as call to
!   this routine occurrs on the first output timestep.
!
              IF (FILETYPE == 'p') THEN
                FILETYPE_2 = LETTER_3
                IF (STYLE /= 'C') THEN
                  IF (ISUBMODEL == atmos_sm) THEN
                    IF (STEPim(a_im) == 0) THEN
                     STEPS = REINIT_STEPS
                    ELSE
                     STEPS = STEPim(a_im) + REINIT_STEPS - 1
                    ENDIF
! DEPENDS ON: stp2time
                    CALL STP2TIME(STEPS,                                &
     &                            A_STEPS_PER_HR*24,86400,              &
     &                            END_DAYS,END_SECONDS)
                  ELSE IF (ISUBMODEL == ocean_sm) THEN
                    IF (STEPim(o_im) == 0) THEN
                     STEPS = REINIT_STEPS
                    ELSE
                     STEPS = STEPim(o_im) + REINIT_STEPS - 1
                    ENDIF
! DEPENDS ON: stp2time
                    CALL STP2TIME(STEPS,                                &
     &                STEPS_PER_PERIODim(o_im),SECS_PER_PERIODim(o_im), &
     &                            END_DAYS,END_SECONDS)
                  ELSE IF (ISUBMODEL == wave_sm) THEN
                    IF (STEPim(w_im) == 0) THEN
                     STEPS = REINIT_STEPS
                    ELSE
                     STEPS = STEPim(w_im) + REINIT_STEPS - 1
                    ENDIF
! DEPENDS ON: stp2time
                    CALL STP2TIME(STEPS,                                &
     &                STEPS_PER_PERIODim(w_im),SECS_PER_PERIODim(w_im), &
     &                            END_DAYS,END_SECONDS)
                  ENDIF
! DEPENDS ON: sec2time
                  CALL SEC2TIME(END_DAYS,END_SECONDS,                   &
     &                          BASIS_TIME_DAYS,BASIS_TIME_SECS,        &
     &                          YYYY,MM,DD,HH,IMIN,ISEC,DAYNO,LCAL360)
                ENDIF
              ENDIF
            ENDIF

          ELSE         !  MEANS
!
!L determine if special mean period
            IF (ISUBMODEL == ATMOS_SM) THEN
              HOURS=DUMPFREQim(a_im)/A_STEPS_PER_HR
              DO I=1,MEANLEV
                HOURS=HOURS*MEANFREQim(I,a_im) !hours per meaning period
              ENDDO
            ELSE IF (ISUBMODEL == OCEAN_SM) THEN
              HOURS=DUMPFREQim(o_im)*SECS_PER_PERIODim(o_im)            &
     &                        /(3600*STEPS_PER_PERIODim(o_im))
              DO I=1,MEANLEV
                HOURS=HOURS*MEANFREQim(I,o_im) !hours per meaning period
              ENDDO
            ELSE IF (ISUBMODEL == WAVE_SM) THEN
              HOURS=DUMPFREQim(w_im)*SECS_PER_PERIODim(w_im)            &
     &                        /(3600*STEPS_PER_PERIODim(w_im))
              DO I=1,MEANLEV
                HOURS=HOURS*MEANFREQim(I,w_im) !hours per meaning period
              ENDDO
            ENDIF
            IF (MOD(HOURS,24) == 0) THEN
              DAYS=HOURS/24
! DEPENDS ON: day2char
              CALL DAY2CHAR(DAYS,FILETYPE_2)
              IF (FILETYPE_2 == '0') THEN
                FILETYPE_2=MEAN_PERIOD(MEANLEV)
              ELSE if (FILETYPE_2 == 'm'.and.DD == 1                    &
     &              .AND.TIME_CONVENTION /= 'Absolute_short    ') THEN
                STYLE='C'      ! period starts at beginning of a month
                IF (MM == 1) THEN ! correct year if month december
                  YYYY=YYYY-1
                ENDIF
              ELSE if (FILETYPE_2 == 's'.and.DD == 1                    &
     &               .AND.TIME_CONVENTION /= 'Absolute_short    ') THEN
                STYLE='C'      ! period starts at beginning of aseason
              ENDIF
            ELSE
              FILETYPE_2=MEAN_PERIOD(MEANLEV)
            ENDIF
          ENDIF
!
          Y_UNITS = CHAR_ID(MOD(YYYY,10)+1)
          M = CHAR_ID(MM+1)
          D = CHAR_ID(DD+1)
          H = CHAR_ID(HH+1)
!
        ENDIF
!
      ELSE
! partial sum files - no date time information required
        SEPARATOR='_'
        FILETYPE_2 = MEAN_PERIOD(MEANLEV)
      ENDIF
!
!L----------------------------------------------------------------------
!L 3.1 Construct filename from the various components
!L
      FILENAME="              "
      IF (MODEL_STATUS == 'Operational'.OR.MODEL_STATUS == 'SCS') THEN
        FILENAME(1:2)  =QW
        FILENAME(3:6)  =EXPT_ID
        FILENAME(7:7)  =SUBMODEL_ID
        FILENAME(8:8)  =SEPARATOR
        FILENAME(9:9)  =FILETYPE
        FILENAME(10:10)=FILETYPE_2
        FILENAME(11:11)=HUNDREDS
        FILENAME(12:12)=TENS
        FILENAME(13:13)=UNITS
        IF (TIME_CONVENTION == 'Sub-hourly      ') THEN
          FILENAME(14:14)=DECI
        END IF
      ELSE
        FILENAME(1:4)  =EXPT_ID
        FILENAME(5:5)  =JOB_ID
        FILENAME(6:6)  =SUBMODEL_ID
        FILENAME(7:7)  =SEPARATOR
        FILENAME(8:8)  =FILETYPE
        FILENAME(9:9)  =FILETYPE_2
        IF (FILETYPE == 's') THEN
          IF (TOGGLE == 1) THEN
            FILENAME(10:10)='a'
          ELSE
            FILENAME(10:10)='b'
          ENDIF
        ELSE IF (STYLE == 'A') THEN
          FILENAME(10:10)=HUNDREDS
          FILENAME(11:11)=TENS
          FILENAME(12:12)=UNITS
          IF (TIME_CONVENTION == 'Sub-hourly      ') THEN
            FILENAME(13:13)=DECI
          END IF
        ELSE IF (STYLE == 'B') THEN
          IF (TIME_CONVENTION == 'Absolute_standard') THEN
!
! decades meansured relative to 1800
            Y_TENS  = CHAR_ID(MOD(YYYY/10,36)+1)
            FILENAME(10:10)=Y_TENS
            FILENAME(11:11)=Y_UNITS
            FILENAME(12:12)=M
            FILENAME(13:13)=D
            FILENAME(14:14)=H
          ELSE IF (TIME_CONVENTION == 'Absolute_long    ') THEN

! centuries  measured from 0 ie 1992  as j92, with wraparound at 3600
            Y_HUNDS = CHAR_ID(MOD(YYYY/100,36)+1)
            Y_TENS  = CHAR_ID((YYYY-(YYYY/100)*100)/10+1)
            FILENAME(10:10)=Y_HUNDS
            FILENAME(11:11)=Y_TENS
            FILENAME(12:12)=Y_UNITS
            FILENAME(13:13)=M
            FILENAME(14:14)=D
          ELSE IF (TIME_CONVENTION == 'Absolute_short   ') THEN
            FILENAME(10:10)=Y_UNITS
    1       FORMAT (I3)
    2       FORMAT ('0',I2)
    3       FORMAT ('00',I1)
            IF (DAYNO <  100) THEN
              IF (DAYNO <  10) THEN
                WRITE(CDAYNO,3) DAYNO
              ELSE
                WRITE(CDAYNO,2) DAYNO
              ENDIF
            ELSE
              WRITE(CDAYNO,1) DAYNO
            ENDIF
            FILENAME(11:13)= CDAYNO
            FILENAME(14:14)=H
          ENDIF
        ELSE  ! style C - Character date
          IF (TIME_CONVENTION == 'Absolute_standard') THEN
!
! decades meansured relative to 1800
            Y_TENS  = CHAR_ID(MOD(YYYY/10,36)+1)
            FILENAME(10:10)=Y_TENS
            FILENAME(11:11)=Y_UNITS
            IF (MEANLEV == 0) THEN
              FILENAME(12:14) = MONTH_3CHAR(MM)
            ELSE ! means date routine called is at beginning of next
!                   period
             MON=MM-1
             IF (MON == 0) THEN
               MON = 12
             ENDIF
             IF (FILETYPE_2 == 'm') THEN
               FILENAME(12:14) = MONTH_3CHAR(MON)
             ELSE
               FILENAME(12:14) = SEASON_3CHAR(MON)
             ENDIF
            ENDIF
          ELSE IF (TIME_CONVENTION == 'Absolute_long    ') THEN

! centuries  measured from 0 ie 1992  as j92, with wraparound at 3600
            Y_HUNDS = CHAR_ID(MOD(YYYY/100,36)+1)
            Y_TENS  = CHAR_ID((YYYY-(YYYY/100)*100)/10+1)
            FILENAME(10:10)=Y_HUNDS
            FILENAME(11:11)=Y_TENS
            FILENAME(12:12)=Y_UNITS
            IF (MEANLEV == 0) THEN
              FILENAME(13:14) = MONTH_2CHAR(MM)
            ELSE ! means date routine called is at beginning of next
!                   period
             MON=MM-1
             IF (MON == 0) THEN
               MON = 12
             ENDIF
             IF (FILETYPE_2 == 'm') THEN
               FILENAME(13:14) = MONTH_2CHAR(MON)
             ELSE
               FILENAME(13:14) = SEASON_2CHAR(MON)
             ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
!
 999  CONTINUE
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE GET_NAME
