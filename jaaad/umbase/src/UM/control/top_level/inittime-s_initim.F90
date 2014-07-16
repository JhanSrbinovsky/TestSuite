#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  INITTIME
!!!
!!!  Purpose: Initialises the model time relative to the calender zero
!!!           time. The model basis time (MBT) is converted to a time
!!!           in seconds since T=0 with respect to the calender. The
!!!           model initialisation time in days and seconds is then
!!!           calculated relative to this.
!!!
!!!  INPUTS - Year
!!!         - Month
!!!         - Day
!!!         - Hour
!!!         - Minute
!!!         - Second
!!!         - Flag for 360 day year calender
!!!
!!!  OUTPUTS - Daynumber in year
!!!          - Time within daynumber
!!!
!!!
!!! Phil Hopwood <- programmer of some or all of previous code changes
!!!
!!!
!!!
!!! Documentation: SSFM Project Documentation - Implementation of UM
!!!                date structure within the SCM
!!!
!!!-------------------------------------------------------------------

!!  Arguments --------------------------------------------------------
!
      Subroutine inittime(iyear,imonth,iday,ihour,imin,isec,            &
     &                    dayno,daytime,lcal360)
!
      Implicit none
!
      external time2sec
!
      Integer                                                           &
     & iyear                                                            &
                                   ! IN Year
     &,imonth                                                           &
                                   ! IN Month
     &,iday                                                             &
                                   ! IN Day
     &,ihour                                                            &
                                   ! IN Hour
     &,imin                                                             &
                                   ! IN Minute
     &,isec                        ! IN Second
      Logical                                                           &
     & lcal360                     ! IN Flag for 360 year calender
      Integer                                                           &
     & dayno                                                            &
                                   ! OUT Daynumber in year
     &,daytime                     ! OUT Time in daynumber
      Integer                                                           &
     & basis_time_days                                                  &
                                   ! LOCAL whole days to basis time
!                                    from start of calender
     &,basis_time_secs             ! LOCAL secs in day at basis time
!
! Calculate Model Basis Time (i.e. begining of current year)
! relative to calender zero
!
! DEPENDS ON: time2sec
      Call time2sec(iyear-1,12,31,0,0,0,                                &
     &  0,0,basis_time_days,basis_time_secs,lcal360)
!
! Calculate daynumber and initial time relative to Model Basis Time
!
! DEPENDS ON: time2sec
      Call time2sec(iyear,imonth,iday,ihour,imin,isec,                  &
     &  basis_time_days,basis_time_secs,                                &
     &  dayno,daytime,lcal360)
                                !
      Return
      END SUBROUTINE inittime
!
#endif
