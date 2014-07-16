#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE BLRAND(pos,bl,orog,z_top_of_model,                     &
     &  first_constant_r_rho_level,nfill,seed2)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods :   To reposition cells in B.L.
!-
!-   Inputs  : POS,BL
!-   Outputs : POS
!-   Controls:
!-
!
! History:
! Version   Date                    Comment
!  5.4    20/11/02  Created.  W.J. Collins
!  5.5    21/01/04  New CONTAINed function RTOETAI added, as functions
!                   RTOETA and GETMETPOINT cannot be inlined together.
!                   K. Ketelsen.
!  6.1    23/09/04  Now uses integer random number seeds. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!-
!VVV  V5.4  BLRAND 20/XI/02
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: first_constant_r_rho_level
      INTEGER, INTENT(IN) :: nfill

      INTEGER, DIMENSION(ransize), INTENT(INOUT) :: seed2
      REAL,                           INTENT(IN) :: z_top_of_model
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: bl
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: orog
      REAL, DIMENSION(4,nclprc),      INTENT(INOUT) :: pos

      INTEGER :: j

      REAL    :: zorog
      REAL    :: blplus
      REAL    :: zbl
      REAL, DIMENSION(nfill) :: xran

      DO j = 1, nfill          ! Generate a sequence of random numbers
        xran(j) = FRANV(seed2)
      END DO

      DO j=1,nfill

! BL now in height
! DEPENDS ON: getmetpoint
        zorog = GETMETPOINT(pos(:,j),orog,.TRUE.)
! DEPENDS ON: getmetpoint
        zbl = GETMETPOINT(pos(:,j),bl,.TRUE.)

        IF (pos(4,j)-earth_radius-zorog <= zbl) THEN   ! In BL
! random reassignment in bl+extra bit:
          blplus = zbl + z_blplus
          pos(4,j) = xran(j) * blplus + earth_radius + zorog
        END IF

! Retrieve cells above domain.
        IF (pos(4,j) - earth_radius - z_top_of_model > -0.1)            &
     &    pos(4,j) = earth_radius + z_top_of_model - 0.1

! Convert to eta
        pos(3,j) = RTOETAI(pos(4,j),zorog,z_top_of_model,               &
     &    first_constant_r_rho_level)
      END DO

!kk   For some reason GETMETPOINT and RTOETA cannot be expanded
!kk   together. RTOETAI with CONTAINS included here

      CONTAINS

      REAL FUNCTION RTOETAI(r,zorog,z_top_of_model,                     &
     &  first_constant_r_rho_level)
! ----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
! ----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: first_constant_r_rho_level

      REAL, INTENT(IN) :: z_top_of_model
      REAL, INTENT(IN) :: zorog
      REAL, INTENT(IN) :: r

      REAL :: a
      REAL :: b
      REAL :: c
      REAL :: r_a
      REAL :: er_fcr

      r_a = r - earth_radius
      er_fcr = eta_rho(first_constant_r_rho_level)

      RTOETAI = r_a / z_top_of_model

      IF (zorog > 1.0e-05 .AND. RTOETAI < er_fcr) THEN
        a = zorog / (er_fcr ** 2)
        b = z_top_of_model - (2.0*zorog / er_fcr)
        c = zorog - r_a
        RTOETAI = (-b + SQRT(b**2 - 4.0*a*c)) / (2.0*a)
      END IF

      END FUNCTION RTOETAI

      END SUBROUTINE BLRAND
#endif
