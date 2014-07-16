
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE pole_bearing(row_length, rows, lat_rot_NP, &
    long_rot_NP, true_longitude, f3_at_u, two_Omega,  &
    bear_rot_NP)

  IMPLICIT NONE

! Description:
!   Calculates bearing of Grid North (from True North) for gridpoints
!   in a Local Area Model rotated grid.
!
! Method:
!   Spherical trig.
!
! Current Code Owner: James Manners
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 6.2       13/06/05  Original code.  James Manners
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):

! Subroutine arguments

  INTEGER, INTENT(IN) :: &
       row_length, rows         ! grid size
  REAL, INTENT(IN) ::  &
       lat_rot_NP,     &        ! Real latitude and longitude
       long_rot_NP,    &        !  of 'pseudo' N pole in radians.
       two_Omega

  REAL, DIMENSION(row_length,rows), INTENT(IN) :: &
       true_longitude, &        ! true longitude of point (rads)
       f3_at_u

  REAL, DIMENSION(row_length,rows), INTENT(OUT) :: &
       bear_rot_NP              ! Bearing of 'pseudo' N pole (rads)

! Local constants

! Local variables
  REAL, DIMENSION(row_length,rows) :: &
       sin_true_latitude, true_lat, long_diff

! Function definitions

!- End of header

! Calculate true latitude at each grid-point.

  sin_true_latitude = f3_at_u / two_Omega
  true_lat = asin(sin_true_latitude)

! Calculate bearing of Grid North

  long_diff = true_longitude - long_rot_NP

  bear_rot_NP = atan2( -cos(lat_rot_NP)*sin(long_diff),  &
        cos(true_lat)*sin(lat_rot_NP) -                  &
        sin_true_latitude*cos(lat_rot_NP)*cos(long_diff) )

END SUBROUTINE pole_bearing
