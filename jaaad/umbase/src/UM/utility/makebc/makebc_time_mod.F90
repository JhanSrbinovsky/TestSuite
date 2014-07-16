#if defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! *********************************************************************
!+ Module makebc_time_mod : Contains a type which is used by a number of
!                           routines to store time values.

MODULE makebc_time_mod
 
TYPE time
  INTEGER :: year
  INTEGER :: month
  INTEGER :: day
  INTEGER :: hour
  INTEGER :: min
  INTEGER :: sec
  INTEGER :: day_no
END TYPE time

END MODULE makebc_time_mod
#endif
