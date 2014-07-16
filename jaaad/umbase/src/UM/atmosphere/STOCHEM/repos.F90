#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE REPOS(pos,bl,orog,z_top_of_model,                      &
     &                 first_constant_r_rho_level)
!-----------------------------------------------------------------------
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.1    03/10/96  Created.  C.E. Johnson   From old NEWPOS
!  5.3    26/06/01  New Dynamics version. C.E. Johnson.
!  5.5    22/01/04  Function RTOETA inlined by hand. K. Ketelsen
!  6.1    22/08/04  Minor tidying of code. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: first_constant_r_rho_level

      REAL,                           INTENT(IN) :: z_top_of_model
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: bl
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: orog
      REAL, DIMENSION(4),          INTENT(INOUT) :: pos

      REAL :: zorog
      REAL :: zbl

      REAL :: a
      REAL :: b
      REAL :: c
      REAL :: r_a

! BL now in height
! DEPENDS ON: getmetpoint
      zorog = GETMETPOINT(pos,orog,.TRUE.)
! DEPENDS ON: getmetpoint
      zbl = GETMETPOINT(pos,bl,.TRUE.)

! If cell is below surface reposition in mid BL.
      IF (pos(4)-earth_radius-zorog < 0.1)                              &
     &  pos(4)=earth_radius+zorog+zbl/2.0

! Retrieve cells above domain.
      IF (pos(4)-earth_radius-z_top_of_model > -0.1)                    &
     &   pos(4) = earth_radius+z_top_of_model-0.1

! Convert to eta
      r_a = pos(4) - earth_radius
      pos(3) = r_a / z_top_of_model

      IF (zorog>1.0E-05.AND.pos(3)<eta_rho(first_constant_r_rho_level)) &
     &  THEN
        a = zorog / (eta_rho(first_constant_r_rho_level))**2
        b = z_top_of_model - (2.0*zorog/                                &
     &    eta_rho(first_constant_r_rho_level))
        c = zorog - r_a
        pos(3) = (-b+SQRT(b**2 - 4.0*a*c))/(2.0*a)
      END IF

      END SUBROUTINE REPOS
#endif
