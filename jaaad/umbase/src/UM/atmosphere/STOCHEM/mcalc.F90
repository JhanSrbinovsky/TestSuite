#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE MCALC(mass,total,xx,clist,ipos,pos,tropz,nchem,nfill)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : CALCULATES MASS (IN NO. OF MOLECULES) OF
!-                         SPECIES IN EACH GRID VOLUME
!-
!-   Inputs  : XX,IPOS,CLIST,NCHEM,TROPZ
!-   Outputs : MASS
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.4    10/05/94  Created.  W.J. Collins
!  5.1    03/12/01  Updated to new dynamics. W.J. Collins.
!  5.5    23/02/04  Vectorised code. K. Ketelsen
!  6.1    21/10/04  No change.
!VVV  V2.2  MCALC 20/X/99
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      USE IN_STOCHEM_OUT
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER,                                   INTENT(IN) :: nfill
      INTEGER,                                   INTENT(IN) :: nchem
      INTEGER, DIMENSION(3,nclprc),              INTENT(IN) :: ipos
      INTEGER, DIMENSION(numchem),               INTENT(IN) :: clist

      REAL, DIMENSION(nc,nclprc),                INTENT(IN) :: xx
      REAL, DIMENSION(4,nclprc),                 INTENT(IN) :: pos
      REAL, DIMENSION(nlonpe,nlatpe),            INTENT(IN) :: tropz
      REAL, DIMENSION(numchem,nlnpe,nlpe,nlev), INTENT(OUT) :: mass
      REAL, DIMENSION(2,nc),                    INTENT(OUT) :: total

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: l
      INTEGER :: n
      INTEGER :: im
      INTEGER :: jm
      REAL, DIMENSION(nfill) :: meantrop
!
      mass = 0.0
      total = 0.0
!
! CDIR nodep
      DO n=1,nfill
        i=ipos(1,n)-lndat+1
        j=ipos(2,n)-ltdat+1
        k=ipos(3,n)
        mass(1:nchem,i,j,k) = mass(1:nchem,i,j,k) + xx(clist(1:nchem),n)
      END DO

! Total concentration in atmosphere
      DO i = 1, SIZE(total,2)
        DO n = 1, nfill
          total(1,i) = total(1,i) + xx(i,n)    ! all cells
        END DO
      END DO

! Find mean tropopause height by linear interpolation
      DO n = 1, nfill
        l = 1 + pos(2,n) / 10
        meantrop(n) = mean_tropz(l) + (pos(2,n)-10.0*(l-1)) *           &
     &    (mean_tropz(l+1) - mean_tropz(l)) / 10.0
      END DO

! Total below tropopause (doesn't account for orography)
      DO n = 1, nfill
        IF (pos(4,n) < meantrop(n) + earth_radius)                      &
     &    total(2,:) = total(2,:) + xx(:,n)
      END DO

      END SUBROUTINE MCALC
#endif
