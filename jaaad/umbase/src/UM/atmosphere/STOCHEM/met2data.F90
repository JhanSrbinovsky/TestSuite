#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE MET2DATA(yy,x,nlevelsin,nlevelsout,lvert)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Convert field from met grid to data grid
!-
!-   Inputs  : x,nlevelsin,nlevelsout,lvert(optional)
!-   Outputs : yy
!-   Controls:
!-
!
! Current Owner of Code: M.G. Sanderson
!
! History:
! Version   Date                    Comment
!  4.2    08/08/96  Created. W.J. Collins
!  5.3    29/11/01  Added option of vertical interpolation. W.J. Collins
!  5.5    20/02/04  Changed /2.0 to *0.5 to cure bit comparability error
!                   Plus other changes for optimisation. K. Ketelsen.
!  6.1    06/08/04  No change
!
!-
!VVV  V2.2  MET2DATA 20/X/99
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTERFACE
! DEPENDS ON: st_height
        INTEGER FUNCTION ST_HEIGHT(pos,eta_array)
          REAL, INTENT(IN) :: pos
          CHARACTER(*), INTENT(IN) :: eta_array
        END FUNCTION ST_HEIGHT
      END INTERFACE

      INTEGER, INTENT(IN) :: nlevelsin
      INTEGER, INTENT(IN) :: nlevelsout
      REAL, DIMENSION(nlonpe,nlatpe,nlevelsin), INTENT(IN) :: x
      REAL, DIMENSION(nlnpe,nlpe,nlevelsout),  INTENT(OUT) :: yy
      LOGICAL, OPTIONAL, INTENT(IN) :: lvert ! True for rho levels

      INTEGER :: k
      INTEGER :: kk
      INTEGER :: ii
      INTEGER :: jj
      INTEGER :: i
      INTEGER :: j
      INTEGER :: i2

      REAL :: aistep
      REAL :: ajstep
      REAL :: fi
      REAL :: fj
      REAL :: fk
      REAL :: area
      REAL :: sum
      REAL :: elo
      REAL :: ehi
      REAL, DIMENSION(nlev) :: delta_eta_stoch
      REAL, DIMENSION(nlnpe,nlpe,nlevelsin) :: y

      aistep = REAL(nmetlong) / REAL(nlong)
      ajstep = REAL((nmetlat-1)) / REAL(mnlat)
      area = aistep * ajstep

      y = 0.0
      yy = 0.0

      delta_eta_stoch = eta_stochem(1:nlev) - eta_stochem(0:nlev-1)
      delta_eta_stoch = 1.0 / delta_eta_stoch

! Take each Met Grid point (noting that they are at the corners of
! the grid squares not the centres)
      DO j=1, nlatpe
        jj = FLOOR((REAL(j)+lobound-1-.5)/ajstep) ! index on data grid
        IF (jj < 0 .OR. jj > mnlat) CYCLE
        DO i=1, nlonpe
          ii = FLOOR((REAL(i)+lnbound-1-.5)/aistep) ! index on data grid
          IF (ii < 0 .OR. ii > nlong-1) CYCLE
          i2 = ii + 1
! split each grid point into 4, with fractions governed by distance to
! edge of data grid. if more than a grid length away set fractions to 0.
          fi = MAX(ii*aistep+1.5-(i+lnbound-1),0.0)
          fj = MAX(jj*ajstep+1.5-(j+lobound-1),0.0)
          IF (ii == 0) ii = nlong
! add each fraction to the appropriate data grid square and divide by
! the relative areas of the grid squares (assuming rectangular grids)
! Only consider squares where PROCMAP=mype
          IF (jj > 0) THEN    ! not at north pole
            IF (procmap(ii,jj) == mype)                                 &
     &        y(ii-lndat+1,jj-ltdat+1,1:nlevelsin) =                    &
     &        y(ii-lndat+1,jj-ltdat+1,1:nlevelsin) + x(i,j,:)*fi*fj/area
            IF (procmap(i2,jj) == mype)                                 &
     &        y(i2-lndat+1,jj-ltdat+1,1:nlevelsin) =                    &
     &        y(i2-lndat+1,jj-ltdat+1,1:nlevelsin) +                    &
     &        x(i,j,:) * (1.0-fi) * fj / area
          END IF
          IF (jj < mnlat) THEN   ! not at south pole
            IF (procmap(ii,jj+1) == mype)                               &
     &        y(ii-lndat+1,jj+1-ltdat+1,1:nlevelsin) =                  &
     &        y(ii-lndat+1,jj+1-ltdat+1,1:nlevelsin) +                  &
     &        x(i,j,:) * fi * (1.0-fj) / area
            IF (procmap(i2,jj+1) == mype)                               &
     &        y(i2-lndat+1,jj+1-ltdat+1,1:nlevelsin) =                  &
     &        y(i2-lndat+1,jj+1-ltdat+1,1:nlevelsin) +                  &
     &        x(i,j,:) * (1.0-fi) * (1.0-fj) / area
          END IF
        END DO
      END DO

      IF (nlevelsin > 1 .AND. nlevelsin /= nlevelsout) THEN
        DO k=1, nlevelsin
          IF (lvert) THEN
! Only need k==1 test for Eta_rho levels.
            IF (k == 1) THEN
              elo = 0.0
            ELSE
              elo = (eta_rho(k) + eta_rho(k-1))  * 0.5
            END IF
          ELSE
            elo = (eta_theta(k) + eta_theta(k-1)) * 0.5
          END IF
          IF (k == nmetlev) THEN
            ehi = 1.0
          ELSE
            IF (lvert) THEN
              ehi = (eta_rho(k+1) + eta_rho(k)) * 0.5
            ELSE
              ehi = (eta_theta(k+1) + eta_theta(k)) * 0.5
            END IF
          END IF
          IF (elo > eta_stochem(nlev)) EXIT
! DEPENDS ON: st_height
          kk = ST_HEIGHT(elo,'Eta_stochem')
          fk = MIN((eta_stochem(kk)-elo) / (ehi-elo), 1.0)
          yy(:,:,kk) = yy(:,:,kk) +                                     &
     &      y(:,:,k) * fk * (ehi-elo) * delta_eta_stoch(kk)
          IF (kk < nlev) THEN
            yy(:,:,kk+1) = yy(:,:,kk+1) +                               &
     &       y(:,:,k) * (1.0-fk) * (ehi-elo) * delta_eta_stoch(kk)
          END IF
        END DO
      ELSE
        IF (nlevelsin == 1) THEN
          yy(:,:,1) = y(:,:,1)
        ELSE
          yy(:,:,:) = y(:,:,:)
        END IF
      END IF

      END SUBROUTINE MET2DATA
#endif
