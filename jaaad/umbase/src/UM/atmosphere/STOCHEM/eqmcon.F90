#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EQMCON(asize,tc,ke,kh)
!-   Purpose and Methods : Calculate Henry's law and equilibrium constan
!-
!-   Inputs  : tc
!-   Outputs : ke,kh
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.4    07/01/97  Created.  W.J. Collins
!  5.5    03/11/03  Vectorised code. M.G. Sanderson
!  6.1    20/10/04  No change.
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!-
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: asize
      REAL, DIMENSION(chunk), INTENT(IN) ::    tc
      REAL, DIMENSION(chunk,6), INTENT(OUT) :: ke, kh

      INTEGER :: i                   ! loop count
      INTEGER :: j                   ! loop count
      REAL, DIMENSION(chunk) :: reus
      REAL, PARAMETER :: one_over_tref = 1.0 / 298.0
      REAL :: f

! Set up Arrhenius parameters for calculation of ke and kh
! Activation energies are actually -Ea/R (units: K)
      REAL, DIMENSION(6) :: a_ke =                                      &
     &  (/ 1.8e+05, 1.7e-02, 6.0e-08, 1.8e-05, 4.3e-07, 1.8e-16 /)
      REAL, DIMENSION(6) :: ea_ke =                                     &
     &  (/ -450.0,  2090.0,  1120.0,  -450.0,  -913.0,  -6716.0 /)
      REAL, DIMENSION(6) :: a_kh =                                      &
     &  (/ 1.1e-02, 3.3e+06,7.36e+04, 1.23e+0, 7.5e+01, 3.4e-02 /)
      REAL, DIMENSION(6) :: ea_kh =                                     &
     &  (/  2300.0,  8700.0,  6621.0,  3120.0,  3400.0,  2420.0 /)

! Factor to remove mult/div of constants from inside loop
      f = 1.0e3 * rho_h2o / mh2o

      DO j = 1, asize
        reus(j) = (1.0 / tc(j)) - one_over_tref
      END DO

! Equilibrium constants (ke) for:
! 1. HNO3 = NO3(-) + H(+)
! 2. SO2(aq) = H(+) + HSO3(-)
! 3. HSO3(-) = H(+) + SO3(2-)
! 4. NH3 + H2O = NH4(+) + OH(-)
! 5. CO2 = H(+) + HCO3(-)
! 6. H2O = H(+) + OH(-)

! Henrys law constants kh (mol / (l.atm)) for:
! Aqueous phase equilibria.
! 1. O3
! 2. HNO3
! 3. H2O2
! 4. SO2
! 5. NH3
! 6. CO2

      DO i = 1, 6
        DO j = 1, asize
          ke(j,i) = a_ke(i) * EXP(ea_ke(i) * reus(j))
          kh(j,i) = a_kh(i) * EXP(ea_kh(i) * reus(j))
        END DO
      END DO
      DO j = 1, asize
        ke(j,6) = ke(j,6) * f
      END DO

      END SUBROUTINE EQMCON
#endif
