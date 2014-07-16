#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE STRATMASK(tropz,o3um,orog,o3conc,z_top_of_model,       &
     &  first_constant_r_rho_level)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Interpolate O3 to STOCHEM levels in
!-                         Stratosphere
!-
!-   Inputs  : TROPZ,O3UM,OROG
!-
!-   Outputs : O3CONC
!
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.2    22/02/00  Created. W.J. Collins
!  5.5    12/12/01  New pressure interpolation for ND. W.J. Collins
!  5.5    27/04/04  Extensive changes to be compatible with STRATCALC3,
!                   and vectorisation. C.E. Johnson and M.G. Sanderson
!  6.1    21/10/04  No change.
!
!-
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: first_constant_r_rho_level

      REAL, INTENT(IN) :: z_top_of_model
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(IN) :: o3um
      REAL, DIMENSION(nlonpe,nlatpe),         INTENT(IN) :: tropz
      REAL, DIMENSION(nlonpe,nlatpe),         INTENT(IN) :: orog
      REAL, DIMENSION(nlnpe,nlpe,nlev),       INTENT(INOUT) :: o3conc

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k

      REAL :: rlong
      REAL :: rlat
      REAL :: zorog
      REAL, DIMENSION(4,nlnpe) :: pos2
      REAL, DIMENSION(nlnpe,nlpe,nlev) :: o3_stochem_grid
      LOGICAL, DIMENSION(nlnpe) :: flag

! DEPENDS ON: met2data
      CALL MET2DATA(o3_stochem_grid,o3um,nmetlev,nlev,.FALSE.)

      DO k = nlev, 1, -1
        DO j=1,nlpe
          rlat = (lat(j+ltdat-2) + lat(j+ltdat-1)) / 2.0
          DO i=1,nlnpe
            IF (i+lndat-1 < nlong) THEN
              rlong = (long(i+lndat-1)+long(i+lndat-1+1)) / 2.0
            ELSE
              rlong = (long(i+lndat-1)+                                 &
     &          long(MOD(i+lndat-1,nlong)+1)+360.0) / 2.0
            END IF
            pos2(1,i) = rlong
            pos2(2,i) = rlat
          END DO
          pos2(3,:) = eta_stochem(k)
          pos2(4,:) = 0.0

          DO i=1,nlnpe
! DEPENDS ON: getmetpoint
            zorog = GETMETPOINT(pos2(:,i),orog,.TRUE.)
! flag set to true if location is at or above tropopause
! DEPENDS ON: etator
            flag(i) = (ETATOR(pos2(3,i),zorog,z_top_of_model,           &
     &        first_constant_r_rho_level) -                             &
! DEPENDS ON: getmetpoint
     &        earth_radius-zorog >= GETMETPOINT(pos2(:,i),tropz,.TRUE.))
          END DO
          DO i = 1, nlnpe
            IF (flag(i)) THEN
! Update stratospheric ozone
              o3conc(i,j,k) = o3_stochem_grid(i,j,k)
            END IF
          END DO
        END DO
      END DO

      END SUBROUTINE STRATMASK
#endif
