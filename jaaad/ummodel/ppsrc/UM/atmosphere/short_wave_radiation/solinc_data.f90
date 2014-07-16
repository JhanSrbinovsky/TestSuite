

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Global data module for variables concerned with solar incidence.

MODULE solinc_data

  IMPLICIT NONE
  SAVE

! Description:
!   Global data necessary for calculating the angle of solar incidence
!   on sloping terrain.
!
! Method:
!   Provides global data.
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

  REAL, ALLOCATABLE, DIMENSION(:,:) :: slope_aspect, slope_angle
  REAL, ALLOCATABLE, DIMENSION(:,:) :: sol_bearing, f_orog, orog_corr
  REAL, ALLOCATABLE, DIMENSION(:)   :: lg_orog_corr, lg_f_orog
  LOGICAL :: L_orog = .FALSE.

! slope_aspect: The direction faced by the mean slope - i.e. the
!               bearing of the slope normal projected on the surface
!               measured in radians clockwise from grid north.
!
! slope_angle:  Angle of the mean slope measured in radians from
!               the horizontal.
!
! sol_bearing:  Mean local bearing of the sun over the timestep
!               expressed in radians clockwise from grid north.
!
! orog_corr:    correction factor for the direct solar flux
!               reaching the surface for sloping terrain.
!
! lg_orog_corr: orog_corr gathered at lit points
!
! lg_f_orog:    The extra direct solar flux at the surface due to
!               the orography correction. This is used in the
!               correction to the sw_incs calculation.
!
! f_orog:       lg_f_orog ungathered onto the full grid. This is
!               used to correct net_atm_flux.
!
! L_orog:       model switch for orography scheme
!
!- End of header

END MODULE solinc_data
