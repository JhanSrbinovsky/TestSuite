#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!     Routine to calculate the year in run and actual time and
!     daynumber in year.
!
!     Modification History:
! Version  Date
!  4.5     07/98     SCM integrated as a standard UM configuration
!                    JC Thil.
!  5.3     11/05/01  Updated to include calculation of previous time
!---------------------------------------------------------------------
      Subroutine TIMECALC(year_init, dayno_init, time_init, timestep    &
     &  ,time_string, lcal360, year, day, time_sec, previous_time       &
     &  ,timehr, timemin)
!
      Implicit none

!     Arguments
      Integer                                                           &
     &  year_init                                                       &
                                ! IN year at model start
     &  ,dayno_init             ! IN daynumber in year of start of
                                !  of the model
      Real                                                              &
     &  time_init                                                       &
                                ! IN start time in day of the
                                !  model (seconds).
     &  ,timestep               ! IN no. of secs. between physics
      Logical                                                           &
     &  lcal360                 ! IN true in 360 idealized year

      Real                                                              &
     &  time_sec                ! OUT actual time of day in secs.
                                !    timesteps.
      Integer                                                           &
     &  year                                                            &
                                ! OUT year
     &  ,day                                                            &
                                ! OUT day in year
     &  ,previous_time(7)       ! OUT
      Character*8                                                       &
     &  time_string             ! OUT actual time in day XX..XX..XX
                                !                        hr..mn..se
!     Local Variables
      Integer                                                           &
     &  time,                                                           &
                                ! Seconds elapsed since model ref
                                !  time.
     &  timesec, timemin                                                &
     &  ,timehr                                                         &
     &  ,elapsed_days_prev                                              &
                              ! Days elapsed, end of previous step
     &  ,elapsed_secs_prev                                              &
                              ! Seconds elapsed, end of previous step
     &  ,month, day_m                                                   &
     &  ,elapsed_days                                                   &
                                ! Days elapsed since model reference
                                !  time
     &  ,elapsed_secs                                                   &
                                ! Seconds in day
     &  ,basis_time_days                                                &
                                !
     &  ,basis_time_secs                                                &
                                ! Basis time of the model.
     &  ,timestep_count         ! Counter for timesteps
      Data timestep_count /0/   ! Initialise to 1
!
      time_string = '00.00.00'  ! total secs. up to start
                                !  of this timestep

      time = 24*3600*(dayno_init) + time_init                           &
     &  +    (timestep_count * timestep)
      elapsed_days_prev = int(time/(24*3600))
      elapsed_secs_prev = time - (elapsed_days_prev * (24*3600))

      timestep_count = timestep_count + 1

      time = 24*3600*(dayno_init) + time_init                           &
     &  +    (timestep_count * timestep)
      elapsed_days = int(time/(24*3600))
      elapsed_secs = time - (elapsed_days * (24*3600))

! DEPENDS ON: time2sec
      Call TIME2SEC(year_init-1, 12, 31, 0, 0, 0,                       &
     &  0, 0, basis_time_days, basis_time_secs, lcal360)

! DEPENDS ON: sec2time
      CALL SEC2TIME(ELAPSED_DAYS_PREV,ELAPSED_SECS_PREV,                &
     &              BASIS_TIME_DAYS,BASIS_TIME_SECS,                    &
     &              PREVIOUS_TIME(1),PREVIOUS_TIME(2),PREVIOUS_TIME(3), &
     &              PREVIOUS_TIME(4),PREVIOUS_TIME(5),PREVIOUS_TIME(6), &
     &              PREVIOUS_TIME(7),LCAL360)

! DEPENDS ON: sec2time
      Call SEC2TIME(elapsed_days, elapsed_secs                          &
     &  ,basis_time_days, basis_time_secs                               &
     &  ,year, month, day_m, timehr, timemin, timesec                   &
     &  ,day, lcal360)


!
!     Set up time string for O/P with diagnostics
!
      Write (time_string(1:2), '(i2)') timehr
      Write (time_string(4:5), '(i2)') timemin
      Write (time_string(7:8), '(i2)') timesec
      time_sec = timehr*3600 + timemin*60 + timesec

!
      Return
      END SUBROUTINE TIMECALC
#endif
