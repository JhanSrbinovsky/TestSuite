

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL SUBROUTINE TIMER ------------------------------------------------
!LL
!LL                    Purpose:
!LL Allows the recording of time spent in any section of the program
!LL Two types of timings are supported:
!LL non-inclusive : if a timed section of code (1) contains another
!LL                 timed section of code (2), then the timing for
!LL                 section (1) will not include the time spent in
!LL                 section (2). This is the normal use for the timer
!LL                 routine in the UM up to vn3.4
!LL inclusive     : allows the user to measure the time taken between
!LL                 any two points in the code, irrespective of any
!LL                 other calls to the timer routine within the timed
!LL                 section
!LL
!LL NB: Non-inclusive timers DO INCLUDE any inclusive timer sections
!LL     contained within them. If this section of code should not be
!LL     included, then also time it with a non-inclusive timer
!LL
!LL Timer now also records the time spent in itself
!LL Parameters:
!LL section_name - 20 byte character string containing name of
!LL                timer reference
!LL
!LL action:
!LL  1 -> first call to timer (timer initialisation)
!LL  2 -> last call to timer (prints out the collected data)
!LL  3 -> non-inclusive start timer
!LL  4 -> non-inclusive end timer
!LL  5 -> inclusive start timer
!LL  6 -> inclusive end timer
!LL
!LL Timer should be called with action=1 before the first executable
!LL statement, and with action=2 after the last executable statement.
!LL
!LL   Model               Modification History
!LL  version    Date
!LL   4.1       16/08/94  Based on old UM timer routine
!LL   4.2       08/10/96  Corrected intermediate timer error in
!LL                       elapsed times.
!LL                       Corrected size of message arg. in
!LL                       TIMER_OUTPUT.
!LL                       P.Burton
!LL   4.3       23/01/97  Added better overview for mpp runs
!LL                       Corrected T3E wallclock time calculation
!LL                       P.Burton
!LL   4.5       17/04/98  Added barrier to allow imbalance to be
!LL                       included in the correct timer section
!LL                                                     P.Burton
!LL   4.5       09/07/98  Replaced missing array index for
!LL                       last_ni_wallclock_time_elapsed.
!LL                                                     P. Burton
!LL   5.1       10/04/00  New reconfiguration support. P.Selwood.
!LL   5.2       06/04/99  Moved the initial synchronisation in
!LL                       TIMER_OUTPUT to make sure that the
!LL                       data is not overwritten, if PE 0 is getting
!LL                       from the other PE's.
!LL                                           Author: Bob Carruthers
!LL   5.3       08/06/01  Move get_cpu_time and get_wallclock_time
!LL                       functions into timefn2a.  A van der Wal
!     5.3       26/07/01  Take an internal copy of "action" variable
!                         so as not to modify it (intent in). P.Selwood
!     6.0       12/09/03  Add FLDCALC def. D Robinson
!     6.1       19/10/04  Add MAKEBC def, R. Sempers
!LL
!LL  Author : Paul Burton
!LL

      SUBROUTINE TIMER(section_name,action_arg)


      Use Rcf_parvars_mod, Only :                                       &
     &    mype,                                                         &
     &    nproc,                                                        &
     &    nproc_x,                                                      &
     &    nproc_y,                                                      &
     &    maxproc

      IMPLICIT NONE

! Arguments:
      CHARACTER(LEN=*)    ::  section_name  ! reference name for timed
                                            ! section
      INTEGER, INTENT(IN) ::  action_arg    ! what action to take

! Local variables:
! ni prefix = non-inclusive timings
! in prefix = inclusive timings


      INTEGER action          ! Local copy of action_arg
      INTEGER max_timers      ! maximum number of timings to be handled
         PARAMETER ( max_timers=300 )

      CHARACTER*20 ni_timer_name(max_timers),                           &
                              ! names of timer references
     &             in_timer_name(max_timers)
                              ! names of timer references

      INTEGER ni_number_of_times_timed(max_timers),                     &
                              ! number of times that a section of code
                              ! has been timed
     &        in_number_of_times_timed(max_timers)
                              ! number of times that a section of code
                              ! has been timed

      INTEGER old_timer(max_timers)
                              ! the reference of the timer stopped when
                              ! a new one is started

      REAL ni_cpu_time_elapsed(max_timers),                             &
                              ! total amount of cpu time spent in
                              ! a section of code
     &     ni_wallclock_time_elapsed(max_timers),                       &
                              ! total amount of wallclock time
     &     in_cpu_time_elapsed(max_timers),                             &
                              ! total amount of time spent in a section
                              ! of code
     &     in_wallclock_time_elapsed(max_timers)
                              ! total amount of wallclock time

      REAL ni_cpu_time_started,                                         &
                              ! for non-inclusive timer - cpu time
                              ! of starting
     &     ni_wallclock_time_started,                                   &
                              ! wallclock time of starting
     &     in_cpu_time_started(max_timers),                             &
                              ! for inclusive timer - cpu time
                              !of starting
     &     in_wallclock_time_started(max_timers)
                              ! wallclock time of starting

      INTEGER current_timer   ! for non-inclusive timer - current
                              ! section of code being timed
      INTEGER ni_number_of_timers,                                      &
                              ! number of timers currently known about
     &        in_number_of_timers
                              ! number of timers currently known about

      LOGICAL in_timer_running(max_timers)
                              ! is a particular timer running?

      INTEGER section_ref     ! reference of the current section



      REAL cpu_time_into_timer,                                         &
                              ! cpu time at which timer routine entered
     &     wallclock_time_into_timer ! wallclock "

      LOGICAL timer_on        ! set to FALSE if an error occurs

! Saved variables:
      SAVE ni_timer_name,in_timer_name,                                 &
     &     ni_number_of_times_timed,in_number_of_times_timed,           &
     &     ni_cpu_time_elapsed,in_cpu_time_elapsed,                     &
     &     ni_wallclock_time_elapsed,in_wallclock_time_elapsed,         &
     &     ni_cpu_time_started,in_cpu_time_started,                     &
     &     ni_wallclock_time_started,in_wallclock_time_started,         &
     &     current_timer ,                                              &
     &     ni_number_of_timers, in_number_of_timers,                    &
     &     in_timer_running,                                            &
     &     old_timer, timer_on

! Magic numbers (action types):
      INTEGER first_call_to_timer,                                      &
     &        last_call_to_timer,                                       &
     &        non_incl_start_timer,                                     &
     &        non_incl_end_timer,                                       &
     &        incl_start_timer,                                         &
     &        incl_end_timer,                                           &
     &        intermediate_output

      PARAMETER ( first_call_to_timer = 1,                              &
     &            last_call_to_timer = 2,                               &
     &            non_incl_start_timer = 3,                             &
     &            non_incl_end_timer = 4,                               &
     &            incl_start_timer = 5,                                 &
     &            incl_end_timer = 6,                                   &
     &            intermediate_output = 7 )





      INTEGER info

! Loop counters etc.
      INTEGER I



      EXTERNAL get_cpu_time,get_wallclock_time
      REAL get_cpu_time,get_wallclock_time
! ----------------------------------------------------------------------

      action = action_arg   ! allows local modifications of the argument


      IF (( action  ==  first_call_to_timer) .OR.                       &
     &    ( action  ==  non_incl_start_timer) .OR.                      &
     &    ( action  ==  incl_start_timer)) THEN



        CALL GC_GSYNC(nproc,info)

      ENDIF

      IF (action  >   100) action=action-100
! The following line is useful for general debugging purposes
! It prints out the name of every timed routine on entry and exit
!!!         WRITE(6,*) section_name,' action= ',action

! start up the timer timer

! DEPENDS ON: get_cpu_time
      cpu_time_into_timer = get_cpu_time()
! DEPENDS ON: get_wallclock_time
      wallclock_time_into_timer = get_wallclock_time()
      in_number_of_times_timed(1) = in_number_of_times_timed(1) + 1

! check the length of the section_name

      IF (LEN(section_name)  >   20) THEN
        WRITE(6,*) 'TIMER has detected a non-fatal ERROR'
        WRITE(6,*) 'Section name ',section_name,' is too long.'
        WRITE(6,*) 'Maximum of 20 characters is allowed'
        WRITE(6,*) section_name,' will be truncated to 20 chars.'
      ENDIF


! diagnose what action to take:

      IF (action  ==  first_call_to_timer) THEN

! First call to timer - do initialisation

        DO I=1,max_timers
          ni_timer_name(I)            = '                    '
          in_timer_name(I)            = '                    '

          ni_number_of_times_timed(I) = 0
          in_number_of_times_timed(I) = 0

          ni_cpu_time_elapsed(I)      = 0.
          in_cpu_time_elapsed(I)      = 0.
          ni_wallclock_time_elapsed(I)= 0.
          in_wallclock_time_elapsed(I)= 0.

          in_timer_running(I)=.FALSE.
        ENDDO

        timer_on = .TRUE.

        current_timer = 1
        ni_number_of_timers = 1
        in_number_of_timers = 1
        in_timer_name(1) = 'TIMER'

! and start the timer running

! DEPENDS ON: get_cpu_time
        ni_cpu_time_started = get_cpu_time()
! DEPENDS ON: get_wallclock_time
        ni_wallclock_time_started = get_wallclock_time()
        ni_number_of_times_timed(current_timer) = 1
        old_timer(current_timer) = 0
        ni_timer_name(current_timer) = section_name

! ----------------------------------------------------------------------

      ELSEIF (timer_on .AND.                                            &
     &        ( (action  ==  last_call_to_timer) .OR.                   &
     &          (action  ==  intermediate_output) ) )THEN

! Last call to timer - or intermediate output required, so
! print out table of results

      IF (action  ==  last_call_to_timer) THEN
! the only active timer should be no.1

        IF (current_timer  /=  1) THEN
          WRITE(6,*) 'TIMER has detected an ERROR'
          WRITE(6,*) 'Attempted to print results without switching ',   &
     &             'off all running non-inclusive timers.'
          WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
          timer_on = .FALSE.
          GOTO 9999
        ENDIF

! Make sure there are no inclusive timers still running

        section_ref = 0
        DO I=1,in_number_of_timers
          IF (in_timer_running(I)) section_ref = I
        ENDDO

        IF (section_ref  /= 0) THEN
          WRITE(6,*) 'TIMER has detected an ERROR'
          WRITE(6,*) 'Attempted to print results without switching ',   &
     &             'off all running inclusive timers.'
          WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
          timer_on = .FALSE.
          GOTO 9999
        ENDIF

! Just to make sure that timer isn't called again
        timer_on = .FALSE.

! and switch off the top level non-inclusive timer

        ni_cpu_time_elapsed(current_timer) =                            &
     &      ni_cpu_time_elapsed(current_timer) +                        &
! DEPENDS ON: get_cpu_time
     &      get_cpu_time() - ni_cpu_time_started

        ni_wallclock_time_elapsed(current_timer) =                      &
     &      ni_wallclock_time_elapsed(current_timer) +                  &
! DEPENDS ON: get_wallclock_time
     &      get_wallclock_time() - ni_wallclock_time_started

      ENDIF ! If this is the final call to timer

! DEPENDS ON: timer_output
      CALL TIMER_OUTPUT(                                                &
     &  in_number_of_timers, ni_number_of_timers,                       &
     &  in_cpu_time_elapsed, ni_cpu_time_elapsed,                       &
     &  in_wallclock_time_elapsed, ni_wallclock_time_elapsed,           &
     &  in_number_of_times_timed, ni_number_of_times_timed,             &
     &  in_timer_name, ni_timer_name,                                   &
     &  action,section_name)



! that's it

! ----------------------------------------------------------------------

      ELSEIF (timer_on .AND.                                            &
     &        (action  ==  non_incl_start_timer) ) THEN

! Start a non-inclusive timer running

! Switch off the current timer
      ni_cpu_time_elapsed(current_timer) =                              &
     &    ni_cpu_time_elapsed(current_timer) +                          &
! DEPENDS ON: get_cpu_time
     &    get_cpu_time() - ni_cpu_time_started

      ni_wallclock_time_elapsed(current_timer) =                        &
     &    ni_wallclock_time_elapsed(current_timer) +                    &
! DEPENDS ON: get_wallclock_time
     &    get_wallclock_time() - ni_wallclock_time_started

! See if we're already keeping records for this section

      section_ref = 0
      DO I=1,ni_number_of_timers
        IF (ni_timer_name(I)  ==  section_name) section_ref = I
      ENDDO

! Check to make sure that there is no timer already running for
! this section
! (NB an inclusive timer with the same reference name is allowed
!  to run simultaneously with this one)

      IF (section_ref  ==  current_timer) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Simultaneous non-inclusive timers attempted ',      &
     &             'for section ',section_name
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! calculate the section reference for the new timer

      IF (section_ref  ==  0) THEN
!       this is a new section
        section_ref = ni_number_of_timers+1
        ni_timer_name(section_ref) = section_name
        ni_number_of_timers = section_ref
      ENDIF

! check that max_timers isn't exceeded:
      IF (ni_number_of_timers  >   max_timers) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'More than ',max_timers,' non-inclusive ',           &
     &             'timers is not allowed.'
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! set up old_timer so that when this new timer is stopped, the
! current timer (that we've just stopped) can be restarted

      old_timer(section_ref)=current_timer

! now start up the new timer

      current_timer = section_ref
      ni_number_of_times_timed(current_timer) =                         &
     &  ni_number_of_times_timed(current_timer) + 1
! DEPENDS ON: get_cpu_time
      ni_cpu_time_started = get_cpu_time()
! DEPENDS ON: get_wallclock_time
      ni_wallclock_time_started = get_wallclock_time()

! that's it

! ----------------------------------------------------------------------

      ELSEIF (timer_on .AND.                                            &
     &        (action  ==  non_incl_end_timer) ) THEN

! Stop a non-inclusive timer

! Make sure that we're being asked to end a timer that's actually
! running.

      IF (ni_timer_name(current_timer)  /=  section_name) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Attempted to stop a non-active ',                   &
     &             'non-inclusive timer ',section_name
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! OK - so stop this timer:

      ni_cpu_time_elapsed(current_timer) =                              &
     &    ni_cpu_time_elapsed(current_timer) +                          &
! DEPENDS ON: get_cpu_time
     &    get_cpu_time() - ni_cpu_time_started

      ni_wallclock_time_elapsed(current_timer) =                        &
     &    ni_wallclock_time_elapsed(current_timer) +                    &
! DEPENDS ON: get_wallclock_time
     &    get_wallclock_time() - ni_wallclock_time_started

! and now restart the old timer (ie. the one that was in
! operation at the time this one was started)

      IF (old_timer(current_timer)  ==  0) THEN
! this means I have just stopped the top level timer - there
! are no more to stop. This is an error - I should do this
! by calling the timer with action=2
         WRITE(6,*) 'TIMER has detected an ERROR'
         WRITE(6,*) 'The top-level timer has been stopped'
         WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

      current_timer=old_timer(current_timer)
! DEPENDS ON: get_cpu_time
      ni_cpu_time_started=get_cpu_time()
! DEPENDS ON: get_wallclock_time
      ni_wallclock_time_started=get_wallclock_time()

! ----------------------------------------------------------------------

      ELSEIF (timer_on .AND.                                            &
     &        (action  ==  incl_start_timer) ) THEN

! Start an inclusive timer running

! See if we're already keeping records for this section

      section_ref = 0
      DO I=1,in_number_of_timers
        IF (in_timer_name(I)  ==  section_name) section_ref = I
      ENDDO

! and calculate the section reference

      IF (section_ref  ==  0) THEN
!       this is a new one
        section_ref = in_number_of_timers + 1
        in_timer_name(section_ref) = section_name
        in_number_of_timers = section_ref
      ENDIF

! check that max_timers isn't exceeded:

      IF (in_number_of_timers  >   max_timers) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'More than ',max_timers,' inclusive ',               &
     &             'timers is not allowed.'
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! Check to make sure that there is no timer already running for
! this section
! (NB a non-inclusive timer with the same reference name is allowed
!  to run simultaneously with this one)

      IF (in_timer_running(section_ref)) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Inclusive timer already running for ',              &
     &             section_name
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! so now we can start the timer for this section
      in_number_of_times_timed(section_ref) =                           &
     &  in_number_of_times_timed(section_ref) + 1
      in_timer_running(section_ref) = .TRUE.
! DEPENDS ON: get_cpu_time
      in_cpu_time_started(section_ref) = get_cpu_time()
! DEPENDS ON: get_wallclock_time
      in_wallclock_time_started(section_ref) = get_wallclock_time()

! that's it


! ----------------------------------------------------------------------

      ELSEIF (timer_on .AND.                                            &
     &        (action  ==  incl_end_timer) ) THEN

! Stop an inclusive timer

! Find out what the reference number of this timer is

      section_ref = 0
      DO I=1,in_number_of_timers
        IF (in_timer_name(I)  ==  section_name) section_ref = I
      ENDDO

      IF (section_ref  ==  0) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Attempting to stop a non-existent ',                &
     &             'inclusive timer ',section_name
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! Make sure this timer is actually running at the moment

      IF (.NOT. in_timer_running(section_ref)) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Attempting to stop a non-running ',                 &
     &             'inclusive timer ',section_name
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! now we can stop it
      in_cpu_time_elapsed(section_ref) =                                &
     &   in_cpu_time_elapsed(section_ref) +                             &
! DEPENDS ON: get_cpu_time
     &   get_cpu_time() - in_cpu_time_started(section_ref)

      in_wallclock_time_elapsed(section_ref) =                          &
     &   in_wallclock_time_elapsed(section_ref) +                       &
! DEPENDS ON: get_wallclock_time
     &   get_wallclock_time() - in_wallclock_time_started(section_ref)

      in_timer_running(section_ref) = .FALSE.

! that's it

! ----------------------------------------------------------------------

      ELSEIF (timer_on) THEN

        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Incorrect action= ',action,' supplied by ',         &
     &             'section ',section_name
        WRITE(6,*) 'Non-fatal error - TIMER will continue'

      ENDIF

 9999 CONTINUE

! stop the timer timer
      in_cpu_time_elapsed(1) = in_cpu_time_elapsed(1) +                 &
! DEPENDS ON: get_cpu_time
     &     get_cpu_time() - cpu_time_into_timer
      in_wallclock_time_elapsed(1) = in_wallclock_time_elapsed(1) +     &
! DEPENDS ON: get_wallclock_time
     &     get_wallclock_time() - wallclock_time_into_timer

      RETURN
      END SUBROUTINE TIMER

!*********************************************************************


