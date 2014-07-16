#if defined(A01_3A) || defined(A01_3C) || defined(A01_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Correction for the angle of solar incidence on sloping terrain.
!
! Description:
!   Calculate the orographic correction to be applied to the
!   Direct SW flux at the surface using the angle of solar
!   incidence on sloping terrain.
!   Uses data from global data module 'solinc_data'.
!
! Method:
!   Spherical trigonometry.
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
!

SUBROUTINE solinc(row_length,rows,cos_zenith_angle)

  USE solinc_data, ONLY: sol_bearing, orog_corr, &
                         slope_aspect, slope_angle
  IMPLICIT NONE

! Subroutine arguments

  INTEGER, INTENT(IN) :: &
       row_length, rows         ! grid size
  REAL, DIMENSION(row_length,rows), INTENT(IN) :: &
       cos_zenith_angle         ! Cosine of the solar zenith angle

! Local constants

! Local variables

! Function definitions

!- End of header

! orog_corr is equal to the cosine of the angle between the incoming
! solar insolation and the normal to the mean slope, divided by
! (cosine of the solar zenith angle x cosine of the slope angle).

  WHERE (cos_zenith_angle /= 0.0 .AND. slope_angle /= 0.0)

       orog_corr = 1.0 + TAN( ACOS(cos_zenith_angle) ) * &
          TAN(slope_angle) * COS( sol_bearing - slope_aspect )

       orog_corr = MAX(orog_corr,EPSILON(orog_corr))

  ELSEWHERE
       orog_corr = 1.0
  END WHERE

END SUBROUTINE solinc
#endif
