

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL Subroutine GAS_CALC ----------------------------------------------
!LL
!LL Purpose :
!LL   Calculates the trace gas mixing ratio (or weighting factor for
!LL aerosol forcing fields.  Rates of increase (yearly compound factors)
!LL can be supplied, or spot values (which will be linearly
!LL interpolated) or a mixture of these.  It is designed so it can be
!LL called each time step, but when rates of increase are being used,
!LL values are in fact only updated at New Year.
!LL The rules are:
!LL   If rates exist (i.e. are positive) for the first & current years
!LL then all concentrations are ignored, except for the initial value.
!LL   If there is a positive rate for the current year but not for the
!LL start, the current rate & most recent concentration are used.
!LL   If rates do not exist for the current year then the concentration
!LL is calculated by linear interpolation between the concentrations at
!LL the given years.
!LL   The mixing ratios calculated after the last given year use the
!LL rate for the final given year.
!LL   The 360-day year is assumed.
!LL CARE should be taken if solitary rates are specified, as this can
!LL result in discontinuities in the concentration time profile at
!LL the next given year without a corresponding given rate.
!LL
!LL Authors : Andrew Brady, Tim Johns, William Ingram
!LL
!LL Version for : Cray YMP
!LL
!LL  Model
!LL version  Date
!LL   4.2  19/11/96       Author: William Ingram, reviewer Cath Senior.
!LL   4.4  22/9/97         Correct code for the case when one changes
!LL    from linear interpolation to a rate & then changes the rate.  WJI
! 6.2  21/02/06     Introduced #if defined statements for
!                   version control of radiation code.
!                                  (J.-C.Thelen)
!LL
!LL
!LLEND -----------------------------------------------------------------
!*L Arguments

      SUBROUTINE GAS_CALC(GAS_NOW                                       &
     &                   ,GAS_INDEX_MAX                                 &
     &                   ,GAS_YEAR                                      &
     &                   ,GAS_CONC                                      &
     &                   ,GAS_RATE                                      &
     &                   ,MAX_SCENARIO_PTS                              &
     &                   ,ICODE)

      IMPLICIT NONE

! Maximum size of common block arrays
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
! Submodel parameters for array sizes
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

      REAL          GAS_NOW      !OUT Gas concentration at time step
      INTEGER       GAS_INDEX_MAX                                       &
                                 !IN
     &            , MAX_SCENARIO_PTS ! IN
      INTEGER       GAS_YEAR(MAX_SCENARIO_PTS)   !IN
      REAL          GAS_CONC(MAX_SCENARIO_PTS)   !IN
      REAL          GAS_RATE(MAX_SCENARIO_PTS)   !IN
      INTEGER       ICODE        !OUT Return code: successful=0
!*

! Common blocks

! Model time
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

!     Local variables

      INTEGER       INDEX       ! to subscript gas concs for NOW_TIME
      INTEGER       I           ! Loop over indices
      INTEGER       YEAR_IN_SECS! Year length in seconds
      INTEGER       NOW_TIME_DAY, NOW_TIME_SEC
!                               ! Time now in days/secs from time zero
      INTEGER       GAS_YR_DAY1,  GAS_YR_SEC1
!                               ! Time in days/secs of current GAS_YEAR
      INTEGER       TIME1
!                               ! The same converted to seconds
      INTEGER       GAS_YR_DAY2,  GAS_YR_SEC2
!                               ! Time in days/secs of next GAS_YEAR

!     Check that GASCNST namelist is defined for this year
      IF ( I_YEAR  <   GAS_YEAR(1) ) THEN
        ICODE = 8325
        WRITE (6, *) 'GAS_CALC: no gas data for this year' 
        RETURN
      ENDIF

!     Loop over I to find correct index for current NOW_TIME
      INDEX = 0
      DO I=1, GAS_INDEX_MAX
        IF ( I_YEAR  >=  GAS_YEAR(I) ) INDEX = INDEX+1
      ENDDO

!     Calculate time now in seconds
! DEPENDS ON: time2sec
      CALL TIME2SEC (I_YEAR, I_MONTH, I_DAY, I_HOUR, I_MINUTE, I_SECOND,&
     &              0, 0, NOW_TIME_DAY, NOW_TIME_SEC, .TRUE.)

!     If gas rate at current year is non zero calculate new GAS_NOW
!     by considering compound increases of GAS_RATE(1:INDEX)
      IF ( GAS_RATE(INDEX)  >   0. ) THEN
        YEAR_IN_SECS = 360 * 86400
! DEPENDS ON: time2sec
        CALL TIME2SEC (GAS_YEAR(INDEX), 1, 1, 0, 0, 0,                  &
     &                0, 0, GAS_YR_DAY1, GAS_YR_SEC1, .TRUE.)
        GAS_NOW = GAS_CONC(1)
        DO I=1, INDEX-1
          IF ( GAS_RATE(I)  <   0. ) THEN
             GAS_NOW = GAS_CONC(I+1)
           ELSE
             GAS_NOW = GAS_NOW *                                        &
     &            ( GAS_RATE(I) ** REAL(GAS_YEAR(I+1)-GAS_YEAR(I)) )
          ENDIF
        ENDDO
!       GAS_NOW now holds the concentration in year INDEX - need only
!       update it to the current year.
        GAS_NOW=GAS_NOW*(GAS_RATE(INDEX)**                              &
     &    REAL(((NOW_TIME_DAY-GAS_YR_DAY1)*86400+                       &
     &          NOW_TIME_SEC-GAS_YR_SEC1)/YEAR_IN_SECS))

!     Otherwise calculate by linear interpolation between respective
!     GAS concentrations of given years.
      ELSE
! DEPENDS ON: time2sec
        CALL TIME2SEC (GAS_YEAR(INDEX), 1, 1, 0, 0, 0,                  &
     &                0, 0, GAS_YR_DAY1, GAS_YR_SEC1, .TRUE.)
! DEPENDS ON: time2sec
        CALL TIME2SEC (GAS_YEAR(INDEX+1), 1, 1, 0, 0, 0,                &
     &                0, 0, GAS_YR_DAY2, GAS_YR_SEC2, .TRUE.)
        TIME1   = GAS_YR_DAY1*86400 + GAS_YR_SEC1
        GAS_NOW = GAS_CONC(INDEX) +                                     &
     &          ( GAS_CONC(INDEX+1) - GAS_CONC(INDEX) )                 &
     & * REAL ( NOW_TIME_DAY*86400 + NOW_TIME_SEC - TIME1 )             &
     &      / REAL ( GAS_YR_DAY2*86400 + GAS_YR_SEC2 - TIME1 )
      ENDIF

      RETURN
      END SUBROUTINE GAS_CALC
