
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE aspang_ancil(row_length, rows, land_points, land_sea_mask, &
                        grad_x, grad_y, bear_rot_NP)

  USE solinc_data, ONLY: slope_aspect, slope_angle
  IMPLICIT NONE

! Description:
!   Calculate mean slope aspect and angle in each gridbox and update
!   these variables in the global data module 'solinc_data'
!
! Method:
!   Uses X & Y gradients from ancillary fields.
!
! Current Code Owner: James Manners
!
! Code Description:
!   Language: FORTRAN 95
!
! Declarations:
!
! Global variables (#include statements etc):
!   Data from module 'solinc_data' available by USE association.
!   Include Pi
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------

! Subroutine arguments

  INTEGER, INTENT(IN) :: &
       row_length, rows, &      ! grid size
       land_points              ! number of land points

  LOGICAL, DIMENSION(row_length,rows), INTENT(IN) :: &
       land_sea_mask            ! land-sea mask

  REAL, DIMENSION(land_points), INTENT(IN) :: &
       grad_x,        &         ! orographic X-gradient on land points
       grad_y                   ! orographic Y-gradient on land points

  REAL, DIMENSION(row_length,rows), INTENT(IN) :: &
       bear_rot_NP              ! bearing of 'pseudo' N pole (rads)

! Local variables

  REAL, DIMENSION(land_points) :: &
       work                     ! work array on land points

!- End of header

! Allocate space for arrays in global data module:

  IF (ALLOCATED(slope_angle))  DEALLOCATE(slope_angle)
  IF (ALLOCATED(slope_aspect)) DEALLOCATE(slope_aspect)
  ALLOCATE(slope_angle (row_length,rows))
  ALLOCATE(slope_aspect(row_length,rows))
  slope_angle=0.0
  slope_aspect=0.0

! Find slope angles and aspects from x & y gradients:

  work = ATAN( (grad_x**2 + grad_y**2)**0.5 )
  slope_angle = UNPACK(work, land_sea_mask, slope_angle)

  work = grad_x
  where (grad_x == 0.0) work = EPSILON(grad_x)
  work = Pi - ATAN(grad_y/work) + SIGN(Pi/2.0,work)
  slope_aspect = UNPACK(work, land_sea_mask, slope_aspect)

! Add bearing of 'pseudo' N pole so aspects are relative to
! true North:

  slope_aspect = MODULO(slope_aspect + bear_rot_NP, Pi*2.)

END SUBROUTINE aspang_ancil
