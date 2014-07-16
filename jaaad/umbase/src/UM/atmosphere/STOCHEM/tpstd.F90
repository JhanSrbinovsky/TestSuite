#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE TPSTD(eta,nlevs,tstd,pstd,z_top_of_model,              &
     &  first_constant_r_rho_level)
!
! Calculates standard atmosphere T and P profiles for
! Eta grids
! See UMPD 10, Appendix 2.
!
!+ To calculate standard atmosphere T and P profiles
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  5.5    10/06/03  Created.  C.E. Johnson
!  6.1    21/10/04  Reformatted code. M.G. Sanderson
!
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nlevs
      INTEGER, INTENT(IN) :: first_constant_r_rho_level
      REAL, DIMENSION(nlevs), INTENT(IN) :: eta
      REAL, DIMENSION(nlevs), INTENT(OUT) :: tstd
      REAL, DIMENSION(nlevs), INTENT(OUT) :: pstd
      REAL, INTENT(IN) :: z_top_of_model

      INTEGER :: i
      INTEGER :: k

      REAL :: dz
      REAL :: zin
      REAL :: t

! Z at bottom of layer (m)
      REAL, DIMENSION(6) :: zb = (/0.0E0,1.1E4,2.0E4,3.2E4,4.7E4,5.0E4/)
! P at bottom of layer (Pa)
      REAL, DIMENSION(6) :: pb = (/101325.,22632.,5475.,868.,111.,75./)
! T at bottom of layer (K)
      REAL, DIMENSION(6) :: tb = (/288.15,216.65,216.65,228.65,270.65,  &
     &                             270.65/)
! lapse rate in layer (K/m):
      REAL, DIMENSION(6) :: l = (/-0.0065,0.0,0.001,0.0028,0.0,-0.0028/)

      DO k=1,nlevs
! DEPENDS ON: etator
        zin=ETATOR(eta(k),0.0,z_top_of_model,first_constant_r_rho_level)
        zin = zin - earth_radius
! Find layer
        DO i=1,5
          IF (zin >= zb(i) .AND. zin < zb(i+1)) EXIT
        END DO

        dz = zin - zb(i)
        t = tb(i) + l(i)*dz/2.0
        pstd(k) = pb(i)*exp(-1.0*g*dz*mair/(rmol*t))
        tstd(k) = tb(i) + l(i)*dz
      END DO

      END SUBROUTINE TPSTD
#endif
