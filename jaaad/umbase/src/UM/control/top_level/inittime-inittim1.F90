#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Top Level

SUBROUTINE inittime(                                              &
#include "argduma.h"
  submodel,icode,cmessage)
  !
  !  Routine: INITTIME -------------------------------------------------
  !
  !  Purpose: Initialises the model time relative to the calendar zero
  !           time.  The basis time is converted to a time in seconds
  !           since T=0 with respect to the calendar. If the model basis
  !           time as specified via the history file does not match the
  !           data time in the dump(s), an error is flagged and the
  !           routine exits.
  !           Also sets derived time information from supplied time
  !           information:-
  !           (a) Assimilation start/end timesteps;
  !           (b) Interface field generation start/end steps and length;
  !           (c) Real timestep(s) in seconds.
  !
  !
  !  Programming standard: UM Doc Paper 3, version 8.2 (25/3/2009)
  !
  !  External documentation: On-line UM document C0 - The top-level
  !                          control system
  !

!!$  USE yomhook, ONLY: lhook, dr_hook
!!$  USE parkind1, ONLY: jprb, jpim


  IMPLICIT NONE

  INTEGER      submodel     ! IN  - submodel dump identifier
#include "parparm.h"
#include "typsize.h"
#include "typduma.h"
  INTEGER      icode        ! OUT - Error return code
  CHARACTER*80 cmessage
  !
  !*----------------------------------------------------------------------
  !  Common blocks
  !
#include "cmaxsize.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "chistory.h"
#include "ctime.h"
#include "ctfilt.h"
#include "lbc_coup.h"
#include "clookadd.h"
  !
  !  Local variables
  !
  LOGICAL                                                           &
     l_dt_error              ! Flag for Data Time error
  INTEGER                                                           &
     ical                    ! 1: Gregorian, 2: 360 day calendar
  INTEGER                                                           &
     elapsed_days                                                 &
                                ! Whole days from basis time to VT etc.
     ,    elapsed_secs                                                 &
                                ! Secs-in-day from basis time to VT etc.
     ,    data_minus_basis_days                                        &
                                ! Whole days from basis time to DT
     ,    data_minus_basis_secs                                        &
                                ! Secs-in-day from basis time to DT
     ,    increment_days                                               &
                                ! Days for increment period
     ,    increment_secs                                               &
                                ! Secs for increment period
     ,    final_days                                                   &
                                ! Target in whole days
     ,    final_secs           ! Target in secs
  INTEGER i                                                         &
                                ! Loop counter
     ,       icount               ! job count
  INTEGER :: basis_year
  INTEGER :: basis_month
  INTEGER :: basis_day
  INTEGER :: basis_hour
  INTEGER :: basis_minute
  INTEGER :: basis_second

  INTEGER :: inc_year
  INTEGER :: inc_month
  INTEGER :: inc_day
  INTEGER :: inc_hour
  INTEGER :: inc_minute
  INTEGER :: inc_second

  INTEGER :: target_end_year
  INTEGER :: target_end_month
  INTEGER :: target_end_day
  INTEGER :: target_end_hour
  INTEGER :: target_end_minute
  INTEGER :: target_end_second

  INTEGER :: resubmit_inc_year
  INTEGER :: resubmit_inc_month
  INTEGER :: resubmit_inc_day
  INTEGER :: resubmit_inc_hour
  INTEGER :: resubmit_inc_minute
  INTEGER :: resubmit_inc_second

  INTEGER :: resubmit_target_year
  INTEGER :: resubmit_target_month
  INTEGER :: resubmit_target_day
  INTEGER :: resubmit_target_hour
  INTEGER :: resubmit_target_minute
  INTEGER :: resubmit_target_second

  ! time2sec needs a common reference point of first day of chosen calendar
  ! ie. day number of 1st of January in first year.
  INTEGER,PARAMETER :: day_one_days = 0
  INTEGER,PARAMETER :: day_one_secs = 0

  CHARACTER*5                                                       &
     cdummy               ! used in setting job prefix

  INTEGER                                                           &
     im                                                               &
                                ! internal model id index
     ,ii                                                               &
                                ! internal model loop counter
     ,step                                                             &
     ,steps_per_period                                                 &
     ,secs_per_period                                                  &
     ,a_steps_per_hr                                                   &
     ,group                                                            &
     ,target_end_step                                                  &
     ,h_step                                                           &
     ,h_group

  LOGICAL                                                           &
     newrun                                                           &
                                ! =T: NRUN, ie no internal model steps yet
     ,zerobasis         ! =T: MODEL_BASIS_TIME=0

!!$  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
!!$  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
!!$  REAL(KIND=jprb)               :: zhook_handle
  !
  !----------------------------------------------------------------------
  ! 0. If at start of integration (ie. step numbers both at zero),
  !     set MODEL_BASIS_TIME from validity time of start dump if run is
  !     an assimilation or pseudo-assimilation or if MODEL_BASIS_TIME is
  !     zero (time checking assumed between model time and observations).
  !     Check that timestep definitions are valid.
  !     Set MODEL_DATA_TIME from start dump data time in all cases.
  !
  ! Check for new run (NRUN) or continuation run (CRUN)
!!$  IF (lhook) CALL dr_hook('INITTIME',zhook_in,zhook_handle)

  newrun=.TRUE.
  DO ii=1,n_internal_model
    im=internal_model_list(ii)
    IF(h_stepim(im) /= 0) THEN    ! CRUN, continuation run
      newrun=.FALSE.
    END IF
  END DO ! ii over internal models

  ! Check for MODEL_BASIS_TIME zero
  zerobasis=.TRUE.
  DO i=1,6
    IF(model_basis_time(i) /= 0)THEN
      zerobasis=.FALSE.
    END IF
  END DO

! Compute MODEL_ANALYSIS_HRS from MODEL_ANALYSIS_MINS
  model_analysis_hrs =  REAL(model_analysis_mins)/60.0

  IF (newrun) THEN
    IF ((model_assim_mode == "Atmosphere".OR.                      &
       model_assim_mode == "Coupled   ") .AND.                     &
       (run_assim_mode == "Atmosphere".OR.                         &
       run_assim_mode == "Coupled   "))  THEN
      DO i=1,6
        model_basis_time(i)=a_fixhd(27+i)
      END DO
    ELSE IF ( .NOT.                                                 &
       ((model_assim_mode == "Atmosphere".OR.                      &
       model_assim_mode == "Coupled   ") .AND.                     &
       (run_assim_mode == "Atmosphere".OR.                         &
       run_assim_mode == "Coupled   "))                            &
       .AND.zerobasis) THEN
      DO i=1,6
        model_basis_time(i)=a_fixhd(27+i)
      END DO
    END IF

    DO i=1,6
      model_data_time(i)=a_fixhd(20+i)
      !  Initialise run_resubmit_target to zero
      run_resubmit_target(i)=0
    END DO
    !  Initialise run_resubmit_jobname to current jobname
    run_resubmit_jobname = run_job_name
    !----------------------------------------------------------------------
    ! 1. Check MODEL_BASIS_TIME against dump validity time in header(s)
    !    (first step only - skipped if assimilation or pseudo-assimilation)
    !
    IF ( .NOT.                                                     &
       ((model_assim_mode == "Atmosphere".OR.                      &
       model_assim_mode == "Coupled   ") .AND.                     &
       (run_assim_mode == "Atmosphere".OR.                         &
       run_assim_mode == "Coupled   "))) THEN
      l_dt_error= .NOT.(model_basis_time(1) == a_fixhd(28) .AND.   &
         model_basis_time(2) == a_fixhd(29) .AND.    &
         model_basis_time(3) == a_fixhd(30) .AND.    &
         model_basis_time(4) == a_fixhd(31) .AND.    &
         model_basis_time(5) == a_fixhd(32) .AND.    &
         model_basis_time(6) == a_fixhd(33))
      IF (l_dt_error) THEN
        icode=1
        cmessage="INITTIME: Atmosphere basis time mismatch"
        GO TO 9999
      END IF
    END IF
  ELSE IF (.NOT.newrun .AND. model_basis_time(2) == 0               &
     .AND. model_basis_time(3) == 0 ) THEN

    a_steps_per_hr = 3600*steps_per_periodim(a_im)/                 &
       secs_per_periodim(a_im)
    IF (h_stepim(a_im) >= model_analysis_hrs*a_steps_per_hr) THEN
      !  Convert model_data_time to days/secs
      ! DEPENDS ON: time2sec
      CALL time2sec(model_data_time(1),model_data_time(2),            &
         model_data_time(3),model_data_time(4),                       &
         model_data_time(5),model_data_time(6),                       &
         day_one_days,day_one_secs,elapsed_days,elapsed_secs,lcal360)
      increment_days = 0
      increment_secs = -model_analysis_hrs*3600
      !  Subtract model_analysis_hrs from model_data_time to recover
      !  model_basis_time
      ! DEPENDS ON: time_df
      CALL time_df(elapsed_days,elapsed_secs,increment_days,        &
         increment_secs,final_days,final_secs)
      ! DEPENDS ON: sec2time
      CALL sec2time(final_days,final_secs,day_one_days,day_one_secs,   &
         model_basis_time(1),model_basis_time(2),                      &
         model_basis_time(3),model_basis_time(4),                      &
         model_basis_time(5),model_basis_time(6),                      &
         i_day_number,lcal360)
      WRITE(6,*)' INITTIME; analysis_hrs,incr_secs ',model_analysis_hrs,&
         increment_secs
    ELSE
      DO i = 1,6
        model_basis_time(i) = model_data_time(i)
      END DO
    END IF
    WRITE(6,*)' INITTIME; model_data_time= ',model_data_time
    WRITE(6,*)' INITTIME; model_basis_time= ',model_basis_time
  END IF
  !
  ! 1.2 Check MODEL_DATA_TIME against dump data time in header(s)
  !
  IF (model_status /= 'Operational   ') THEN
    l_dt_error= .NOT.(model_data_time(1) == a_fixhd(21) .AND. &
       model_data_time(2) == a_fixhd(22) .AND.                &
       model_data_time(3) == a_fixhd(23) .AND.                &
       model_data_time(4) == a_fixhd(24) .AND.                &
       model_data_time(5) == a_fixhd(25) .AND.                &
       model_data_time(6) == a_fixhd(26))
    IF (l_dt_error) THEN
      icode=3
      cmessage="INITTIME: Atmosphere data time mismatch"
      GO TO 9999
    END IF
  END IF
  !----------------------------------------------------------------------
  ! 2. Check that model calendar matches the dump header(s)
  !
  IF (lcal360) THEN
    ical=2
  ELSE
    ical=1
  END IF
  IF (a_fixhd(8) /= ical) THEN
    icode=-1
    cmessage="INITTIME: Model calendar doesn't match atmos dump"
    WRITE(6,*)cmessage
    a_fixhd(8)=ical
  END IF
  !
  !  Create copy of logical LCAL360 that can be carried in comdeck ctime
  !   without fear of clashing.
  !
  l_c360dy = lcal360
  !
  !----------------------------------------------------------------------
  ! 3. Initialise model time relative to first day of calendar, using
  !    30 day month calendar for climate models, full calendar including
  !    leap years for other models
  !
  basis_year   = model_basis_time(1)
  basis_month  = model_basis_time(2)
  basis_day    = model_basis_time(3)
  basis_hour   = model_basis_time(4)
  basis_minute = model_basis_time(5)
  basis_second = model_basis_time(6)

  ! DEPENDS ON: time2sec
  CALL time2sec(basis_year,basis_month,basis_day,                       &
     basis_hour,basis_minute,basis_second,                              &
     day_one_days,day_one_secs,basis_time_days,basis_time_secs,lcal360)


  ! 3.1 Set initial time day number in dump header(s) only at start
  !     The i_year outputs from here are not needed later.
  IF(h_stepim(atmos_im) == 0) THEN
    ! DEPENDS ON: sec2time
    CALL sec2time(day_one_days,day_one_secs,basis_time_days,basis_time_secs, &
       i_year,i_month,i_day,i_hour,i_minute,i_second,                        &
       i_day_number,lcal360)
    a_fixhd(27) = i_day_number
  END IF
  !
  !----------------------------------------------------------------------
  ! 4. Initialise incremental step counter(s) to accord with basis time
  !    and model restart time, and flag any difference which might occur
  !    due to altering the timestep partway through an integration.
  !    Set steps-per-group to match timestep, and echo both step and
  !    steps-per-group information in the history file.
  !
  DO ii=1,n_internal_model
    im=internal_model_list(ii)
    secs_per_period =secs_per_periodim(im)
    steps_per_period=steps_per_periodim(im)
    h_step   = h_stepim(im)
    IF(im == atmos_im) THEN
      i_year   = a_fixhd(28)
      i_month  = a_fixhd(29)
      i_day    = a_fixhd(30)
      i_hour   = a_fixhd(31)
      i_minute = a_fixhd(32)
      i_second = a_fixhd(33)
    END IF   ! atmos_im
    !

    ! Calculate elapsed time of run according to dump header.
    ! DEPENDS ON: time2sec
    CALL time2sec(i_year,i_month,i_day,i_hour,i_minute,i_second,   &
       day_one_days,day_one_secs,                            &
       elapsed_days,elapsed_secs,lcal360)
    elapsed_days = elapsed_days - basis_time_days
    elapsed_secs = elapsed_secs - basis_time_secs
    IF (elapsed_secs < 0) THEN
      elapsed_secs = elapsed_secs + 86400
      elapsed_days = elapsed_days - 1
    END IF

    ! DEPENDS ON: time2sec
    CALL time2sec(model_data_time(1),model_data_time(2),   &
       model_data_time(3),model_data_time(4),              &
       model_data_time(5),model_data_time(6),              &
       day_one_days,day_one_secs,                          &
       data_minus_basis_days,data_minus_basis_secs,        &
       lcal360)

    data_minus_basis_days = data_minus_basis_days - basis_time_days
    data_minus_basis_secs = data_minus_basis_secs - basis_time_secs
    IF (data_minus_basis_secs < 0) THEN
      data_minus_basis_secs = data_minus_basis_secs + 86400
      data_minus_basis_days = data_minus_basis_days - 1
    END IF

    data_minus_basis_hrs= data_minus_basis_secs/3600+      &
       24*data_minus_basis_days
    !
    forecast_hrs=       - data_minus_basis_hrs
    !   DATA_MINUS_BASIS_HRS can only be gt 0 for CRUNs when an earlier
    !   run has had MODEL_DATA_TIME updated. For correct values of LBFT
    !   header in CRUNs reset to 0.
    IF(data_minus_basis_hrs >  0) data_minus_basis_hrs=0
    forecast_hrs= elapsed_secs/3600 + 24*elapsed_days                &
       - data_minus_basis_hrs
    IF ((run_assim_mode /= "None      " .OR. l_iau) .AND.            &
       (forecast_hrs >= (model_analysis_hrs-data_minus_basis_hrs)))  &
       forecast_hrs=forecast_hrs-                                    &
       (model_analysis_hrs-data_minus_basis_hrs)
    !
    ! DEPENDS ON: tim2step
    CALL tim2step(elapsed_days,elapsed_secs,                          &
       steps_per_period,secs_per_period,step)
    group= (model_hrs_per_group*3600*steps_per_period)                &
       /secs_per_period
    IF ((elapsed_days /= 0.OR.elapsed_secs /= 0).AND.                 &
       step /= h_step) THEN
      icode=-1
      cmessage="INITTIME: Warning- New STEP doesn't match old value"
      WRITE(6,*)cmessage
      WRITE(6,*)'internal model id',im,' old=',h_step,' New=',step
    END IF

    stepim(im)    = step
    groupim(im)   = group

    h_stepim(im)    = step
    h_groupim(im)   = group
  END DO ! ii over N_INTERNAL_MODEL

  !----------------------------------------------------------------------
  ! 5. Set target end steps from target end time using
  !     relative time convention
  !
  ! Overwrite history resubmit flag with user set switch.
  run_resubmit = control_resubmit
  IF (run_resubmit == "N") THEN

    ! Target_end date = basis date plus run_target_end date
    target_end_year   = basis_year
    target_end_month  = basis_month
    target_end_day    = basis_day
    target_end_hour   = basis_hour
    target_end_minute = basis_minute
    target_end_second = basis_second

    ! DEPENDS ON: add_period_to_date
    CALL add_period_to_date(                                  &
       target_end_year,target_end_month,target_end_day,            &
       target_end_hour,target_end_minute,target_end_second,        &
       run_target_end(1),run_target_end(2),run_target_end(3), &
       run_target_end(4),run_target_end(5),run_target_end(6), &
       lcal360,icode)

    ! Convert the target_end date into the requested end date 
    ! relative to basis date (elapsed date)
    ! DEPENDS ON: time2sec
    CALL time2sec(target_end_year,target_end_month,target_end_day,  &
       target_end_hour,target_end_minute,target_end_second,         &
       day_one_days,day_one_secs,elapsed_days,elapsed_secs,lcal360)

    elapsed_days = elapsed_days - basis_time_days
    elapsed_secs = elapsed_secs - basis_time_secs
    IF (elapsed_secs < 0) THEN
      elapsed_secs = elapsed_secs + 86400
      elapsed_days = elapsed_days - 1
    END IF
    run_resubmit_jobname = ' '
  ELSE
    ! Running with automatic resubmission. Therefore the run length of
    ! this run is the lesser of 
    !        a) the time from now to the run target end.
    !        b) the run resubmit increment

    ! work out next jobname and set present jobname
    !
    run_job_name = run_resubmit_jobname
    READ (run_job_name,'(A5,I3)') cdummy,icount

    icount=MOD(icount+1,1000)

    ! Formatting ensures the numeric part includes leading zeros
    WRITE (run_resubmit_jobname,'(A5,I3.3)') cdummy,icount


    !  run_target_end and run_resubmit_increment are set as a 
    !  period of time, not a date. Need to convert them into 
    !  dates before calculate the actual number of days/seconds 
    !  requested, as. eg. a run period of 1 month is 30 days if
    !  starting in June, but 31 in July.

    !
    !  Work out increment in target date requested by run_resubmit_inc
    !  and convert it into days/seconds since model basis time.
    !
    resubmit_inc_year   = run_resubmit_inc(1)
    resubmit_inc_month  = run_resubmit_inc(2)
    resubmit_inc_day    = run_resubmit_inc(3)
    resubmit_inc_hour   = run_resubmit_inc(4)
    resubmit_inc_minute = run_resubmit_inc(5)
    resubmit_inc_second = run_resubmit_inc(6)

    inc_year   = i_year
    inc_month  = i_month
    inc_day    = i_day
    inc_hour   = i_hour
    inc_minute = i_minute
    inc_second = i_second

    ! DEPENDS ON: add_period_to_date
    CALL add_period_to_date(inc_year,inc_month,inc_day,            &
       inc_hour,inc_minute,inc_second,                             &
       resubmit_inc_year,resubmit_inc_month,resubmit_inc_day,      &
       resubmit_inc_hour,resubmit_inc_minute,resubmit_inc_second,  &
       lcal360,icode)

    ! DEPENDS ON: time2sec
    CALL time2sec(inc_year,inc_month,inc_day,                       &
       inc_hour,inc_minute,inc_second,                              &
       day_one_days,day_one_secs,increment_days,increment_secs,     &
       lcal360)

    ! Calculate time since model basis time
    increment_days = increment_days - basis_time_days
    increment_secs = increment_secs - basis_time_secs
    IF (increment_secs < 0) THEN
      increment_secs = increment_secs + 86400
      increment_days = increment_days - 1
    END IF

    ! Do the same for run_target_end

    target_end_year   = run_target_end(1)
    target_end_month  = run_target_end(2)
    target_end_day    = run_target_end(3)
    target_end_hour   = run_target_end(4)
    target_end_minute = run_target_end(5)
    target_end_second = run_target_end(6)

    ! Add target end period to basis time
    ! DEPENDS ON: add_period_to_date
    CALL add_period_to_date(basis_year,basis_month,basis_day,      &
       basis_hour,basis_minute,basis_second,                       &
       target_end_year,target_end_month,target_end_day,            &
       target_end_hour,target_end_minute,target_end_second,        &
       lcal360,icode)

    ! basis time is now set to target end date. Convert this into
    ! days since Jan 1st year one.
      ! DEPENDS ON: time2sec
    CALL time2sec(basis_year,basis_month,            &
       basis_day,basis_hour,basis_minute,            &
       basis_second,day_one_days,day_one_secs,       &
       final_days,final_secs,lcal360)

    final_days = final_days - basis_time_days
    final_secs = final_secs - basis_time_secs
    IF (final_secs < 0) THEN
      final_secs = final_secs + 86400
      final_days = final_days - 1
    END IF

    ! Choose the shorter of the requested increment and the time to the 
    ! run target end
    IF ((increment_days >  final_days).OR.                                  &
       ((increment_days == final_days).AND.(increment_secs >  final_secs))) &
       THEN
      elapsed_days = final_days
      elapsed_secs = final_secs
    ELSE
      elapsed_days = increment_days
      elapsed_secs = increment_secs
    END IF

  END IF ! End of automatic resubmission section
  IF (elapsed_days <  0) THEN
    icode=1
    cmessage="INITTIME: Negative run length requested"
    GO TO 9999
  END IF
  ! This prevents an infinite resubmission cycle
  IF (MOD(elapsed_secs,secs_per_period/steps_per_period) /= 0)THEN
    icode=1
    cmessage="INITTIME: Run length not integral no. of timesteps. See output"
    WRITE(6,*)'INITTIME: Run length does not divide into timesteps'
    WRITE(6,*)'run length ',elapsed_days,' days ',elapsed_secs,' seconds '
    WRITE(6,*)'TIMESTEP ',secs_per_period/steps_per_period
    WRITE(6,*)'Modify RUN_TARGET_END and resubmit'
    GO TO 9999
  END IF
  !
  DO ii=1,n_internal_model
    im=internal_model_list(ii)
    secs_per_period =secs_per_periodim(im)
    steps_per_period=steps_per_periodim(im)
    ! DEPENDS ON: tim2step
    CALL tim2step(elapsed_days,elapsed_secs,                        &
       steps_per_period,secs_per_period,target_end_step)
    target_end_stepim(im)=target_end_step
  END DO ! ii over N_INTERNAL_MODEL
  !----------------------------------------------------------------------
  ! 7. Set assimilation start timestep, length in steps, and overlap into
  !     forecast from basic control information
  !
  a_steps_per_hr = 3600*steps_per_periodim(a_im)/secs_per_periodim(a_im)
  assim_firststepim(a_im) = ( REAL(a_assim_start_min)/60.0 ) * a_steps_per_hr
  assim_stepsim(a_im) = model_analysis_hrs * a_steps_per_hr - &
     assim_firststepim(a_im)
  assim_extrastepsim(a_im) = &
     ( ( REAL(a_assim_end_min)/60.0 ) -model_analysis_hrs ) * a_steps_per_hr

  !----------------------------------------------------------------------
  ! 7. If running in IAU mode, calculate step on which data time must
  !    be reset. A_STEPS_PER_HR was set in the previous section.

  IF (l_iau .OR. run_assim_mode  ==  "NoIAU     ") THEN

    iau_dtresetstep = model_analysis_hrs * a_steps_per_hr

    IF (h_stepim(a_im) == iau_dtresetstep) THEN
      ! At data time, so reset LBFT to zero for dump fields:
      DO i = 1, a_prog_lookup
        a_lookup(lbft,i)=0
      END DO
    END IF

  END IF
  !----------------------------------------------------------------------
  ! 8. Set LBC_FC_HRS - the number of hours relative to the last data
  !                     time that will apply during this run.

  IF (l_lbc_coup) THEN
    IF (run_assim_mode /= "None      " .OR. l_iau) THEN
      ! Data time reset during run:
      lbc_fc_hrs = -model_analysis_hrs
    ELSE
      ! No data time reset:
      lbc_fc_hrs = forecast_hrs
    END IF
  END IF

  !----------------------------------------------------------------------
  ! 9. Calculate length of timestep in seconds in CTIME
  !
  secs_per_stepim(atmos_sm) = FLOAT(secs_per_periodim(atmos_sm))/ &
     FLOAT(steps_per_periodim(atmos_sm))
  !----------------------------------------------------------------------
  ! 10. Set current time in CTIME according to submodel
  !
  IF (submodel == atmos_sm) THEN
    i_year   = a_fixhd(28)
    i_month  = a_fixhd(29)
    i_day    = a_fixhd(30)
    i_hour   = a_fixhd(31)
    i_minute = a_fixhd(32)
    i_second = a_fixhd(33)
  END IF
  previous_time(1)=i_year
  previous_time(2)=i_month
  previous_time(3)=i_day
  previous_time(4)=i_hour
  previous_time(5)=i_minute
  previous_time(6)=i_second
  previous_time(7)=0        ! Not yet set

9999 CONTINUE

!!$  IF (lhook) CALL dr_hook('INITTIME',zhook_out,zhook_handle)
  RETURN
  !----------------------------------------------------------------------
END SUBROUTINE inittime
#endif
