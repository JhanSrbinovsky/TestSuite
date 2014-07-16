#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: TEMPORAL -------------------------------------------------
!LL
!LL  Purpose: Control routine to handle temporal processing options
!LL           within STASH.  Its input and output arguments look like
!LL           1D arrays (ie. all the data should be in contiguous areas
!LL           of memory).  Lower level service routines are called to
!LL           perform the individual processing options.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   S.Tett
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.1  24/02/93  Change name of variable 'end' to 'last_ts' (ST).
!     4.4  25/11/96  Add processing code option 8 - daily mean
!                    timeseries. R A Stratton.
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D72
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE TEMPORAL(variable,result,size,extra_size,              &
     &  control,control_size,ocean,                                     &
     &  timestep,error,errmssg,start,amdi)
!
      IMPLICIT NONE
!
      INTEGER size                  ! IN  size of arrays
      REAL variable(size)           ! IN  data array
      REAL result(size)             ! OUT output array
      INTEGER extra_size            ! IN size of extra data
      INTEGER control_size          ! IN  size of control
      INTEGER control(control_size) ! IN  control
      INTEGER timestep              ! IN  present value of timestep
      INTEGER error                 ! OUT error code
      CHARACTER*(*) errmssg         ! OUT error message
      REAL amdi                     ! IN  missing data indicator
      LOGICAL ocean                 ! IN  true if ocean diagnostic
      LOGICAL start                 ! OUT true if start timestep
!*----------------------------------------------------------------------
#include "sterr.h"
#include "stparam.h"
!
! Subroutines called
!
      EXTERNAL staccum,stmax,stmin
!
! Local variables
!
      LOGICAL masking        ! indicator for masking (ie. missing data)
      INTEGER proc_code      ! value of processing code
      INTEGER mask_code      ! value of masking code
      REAL divisor           ! divisor for the time mean (1/period)
      INTEGER mod_period     ! timesteps since start modulo period.
      INTEGER start_time     ! value of start time
      INTEGER i              ! loop counter
      LOGICAL last_ts        ! true if end timestep
      INTEGER proc_size      ! size of data to be processed
!L---------------------------------------------------------------------
!L 1. Set processing option code and select appropriate service routine
!L
      proc_size=size-extra_size
      proc_code=control(st_proc_no_code)
!
!  Replace (null processing)
!
      IF (proc_code == st_replace_code) THEN
        DO i=1,size
          result(i)=variable(i)
        ENDDO
        start=(control(st_start_time_code) == timestep)
!
!  Mean/accumulation
!
      ELSEIF (proc_code == st_accum_code.or.                            &
     &        proc_code == st_time_mean_code) THEN
        start_time=control(st_start_time_code)
        IF (control(st_period_code) == st_infinite_time) THEN
          start=(timestep == start_time)
          last_ts=.FALSE.
        ELSE
          mod_period=mod(timestep-start_time,control(st_period_code))
          start=(mod_period == 0)
          last_ts=(mod_period == (control(st_period_code)-              &
     &                        control(st_freq_code)))
        ENDIF
        mask_code=control(st_gridpoint_code)
        mask_code=mod(mask_code,block_size)
        masking=(mask_code /= stash_null_mask_code).or.ocean
        IF (start) THEN      ! first timestep.
          DO i=1,size
            result(i)=variable(i)
          ENDDO
        ELSE
! DEPENDS ON: staccum
          CALL STACCUM(variable,result,proc_size,masking,amdi)
          DO i=proc_size+1,size
            result(i)=variable(i) ! copy over the extra data (if any)
          ENDDO
        ENDIF
!  Normalise at end of mean period
        IF (last_ts.and.proc_code == st_time_mean_code) THEN
          divisor=(float(control(st_freq_code))/                        &
     &             float(control(st_period_code)))
! If field is masked test for MDI, otherwise don't
          IF (masking) THEN
            DO i=1,proc_size
              IF (result(i) /= amdi) THEN
                result(i)=result(i)*divisor
              ENDIF
            ENDDO
          ELSE
            DO i=1,proc_size
              result(i)=result(i)*divisor
            ENDDO
          ENDIF
        ENDIF
!
!  Maximum
!
      ELSEIF (proc_code == st_max_code) THEN
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period == 0)
        IF (start) THEN
          DO i=1,size
            result(i)=variable(i)
          ENDDO
        ELSE
          mask_code=control(st_gridpoint_code)
          mask_code=mod(mask_code,block_size)
          masking=(mask_code /= stash_null_mask_code).or.ocean
! DEPENDS ON: stmax
          CALL STMAX(variable,result,proc_size,masking,amdi)
          DO i=proc_size+1,size
            result(i)=variable(i) ! copy over the extra data (if any)
          ENDDO
        ENDIF
!
!  Minimum
!
      ELSEIF (proc_code == st_min_code) THEN
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period == 0)
        IF (start) THEN
          DO i=1,size
            result(i)=variable(i)
          ENDDO
        ELSE
          mask_code=control(st_gridpoint_code)
          mask_code=mod(mask_code,block_size)
          masking=(mask_code /= stash_null_mask_code).or.ocean
! DEPENDS ON: stmin
          CALL STMIN(variable,result,proc_size,masking,amdi)
          DO i=proc_size+1,size
            result(i)=variable(i) ! copy over the extra data (if any)
          ENDDO
        ENDIF
!
!  Timeseries (append)
!
      ELSEIF (proc_code == st_time_series_code) THEN
        DO i=1,size
! Note that on start timestep this will include the extra data
          result(i)=variable(i)
        ENDDO
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period == 0)
        last_ts=(mod_period == (control(st_period_code)-                &
     &                      control(st_freq_code)))
!
!  Append trajectories
!
      ELSEIF (proc_code == st_append_traj_code) THEN
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period == 0)
        last_ts=(mod_period == (control(st_period_code)-                &
     &                      control(st_freq_code)))
        error=st_not_supported
        write(errmssg,100)' do not support append trajects'
        goto 999
!
!  Timeseries (append) - option 8 daily mean
!
      ELSEIF (proc_code == st_time_series_mean) THEN

        DO i=1,size
! Note that on start timestep this will include the extra data
          result(i)=variable(i)
        ENDDO
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period == 0)
        last_ts=(mod_period == (control(st_period_code)-                &
     &                      control(st_freq_code)))
!
!  Error condition
!
      ELSE
        error=unknown_processing
        write(errmssg,101)' unknown processing code',proc_code
        goto 999
      ENDIF
!
999   CONTINUE   ! jump for errors
!
100   FORMAT('TEMPORAL : >>> FATAL ERROR <<<',a30)
101   FORMAT('TEMPORAL : >>> FATAL ERROR <<<',a30,i5)
!
      RETURN
      END SUBROUTINE TEMPORAL
#endif
