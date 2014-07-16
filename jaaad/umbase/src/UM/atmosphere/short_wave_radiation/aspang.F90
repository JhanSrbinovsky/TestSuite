#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE aspang(row_length, rows, delta_lambda, delta_phi, &
    earth_radius, orog, cos_theta_latitude, bear_rot_NP)

  USE solinc_data, ONLY: slope_aspect, slope_angle
  IMPLICIT NONE

! Description:
!   Calculate mean slope aspect and angle in each gridbox and update
!   these variables in the global data module 'solinc_data'
!
! Method:
!   Method is taken from Gallant & Wilson, Computers & Geosciences
!   Vol. 22, No. 7, pp. 713-722, 1996
!
! Current Code Owner: James Manners
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 6.2       15/02/05  Original code.  James Manners
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):
!   Data from module 'solinc_data' available by USE association.
!   Include Pi
#include "c_pi.h"

! Subroutine arguments

  INTEGER, INTENT(IN) :: &
       row_length, rows         ! grid size
  REAL, INTENT(IN) :: &
       delta_lambda,  &         ! EW (x) grid spacing in radians.
       delta_phi,     &         ! NS (y) grid spacing in radians.
       earth_radius             ! Mean radius of Earth in m.

  REAL, DIMENSION(0:row_length+1,0:rows+1), INTENT(IN) :: &
       orog                     ! Mean gridpoint heights in m.
  REAL, DIMENSION(row_length,rows), INTENT(IN) :: &
       cos_theta_latitude, &    ! Cosine of gridpoint latitude.
       bear_rot_NP              ! Bearing of 'pseudo' N pole (rads)

! Local constants

! Local variables
  REAL, DIMENSION(0:row_length+1,0:rows+1) :: &
       z2,z4,z6,z8              ! Shifted orography grids
  REAL, DIMENSION(row_length,rows) :: &
       zx,            &         ! dz/dx slope in the x direction.
       zy,            &         ! dz/dy slope in the y direction.
       xscale                   ! Width of gridbox (x) in m.
  REAL :: yscale                ! Width of gridbox (y) in m.

! Function definitions

!- End of header

! Allocate space for arrays in global data module:

  IF (ALLOCATED(slope_aspect)) DEALLOCATE(slope_aspect)
  IF (ALLOCATED(slope_angle))  DEALLOCATE(slope_angle)
  ALLOCATE(slope_aspect (row_length,rows), &
           slope_angle  (row_length,rows) )

! Calculate the horizontal scales of the gridboxes (in metres):

  yscale=delta_phi*earth_radius  ! scalar
  xscale=delta_lambda*earth_radius*cos_theta_latitude  ! array
  where (xscale == 0.0) xscale=EPSILON(xscale)

! Calculate x and y slopes, using orography halos at boundaries:

  z2 = EOSHIFT( orog, shift =  1, dim=1 )
  z6 = EOSHIFT( orog, shift = -1, dim=1 )

  z8 = EOSHIFT( orog, shift =  1, dim=2 )
  z4 = EOSHIFT( orog, shift = -1, dim=2 )

  zx=( z2(1:row_length,1:rows) - z6(1:row_length,1:rows) ) / &
                       (2.0*xscale)

  zy=( z8(1:row_length,1:rows) - z4(1:row_length,1:rows) ) / &
                       (2.0*yscale)

  slope_angle  = ATAN( (zx**2 + zy**2)**0.5 )
  where (zx == 0.0) zx=EPSILON(zx)
  slope_aspect = Pi - ATAN(zy/zx) + (Pi/2.0)*( zx/ABS(zx) )

! Add bearing of 'pseudo' N pole so aspects are relative to
! true North.

  slope_aspect = MODULO(slope_aspect + bear_rot_NP, Pi*2.)

END SUBROUTINE aspang
#endif
