#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE DATA2MET2D(x,y,mw)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Convert field from 3-D STOCHEM grid to
!-                         2-D zonal mean on met Eta_theta grid by
!-                         linear interpolation.
!-
!-   Inputs  : X        (Vol. mixing ratio)
!-           : MW       (Molecular Weight of X)
!-   Outputs : Y        (Mass mixing ratio)
!-   Controls:
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  5.5    14/02/03  Created. C.E. Johnson
!  6.1    06/08/04  Reformatted code. M.G. Sanderson
!
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!-----------------------------------------------------------------------

      REAL, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN)  :: x
      REAL, DIMENSION(nlatpe,nmetlev), INTENT(OUT)  :: y
      REAL, INTENT(IN) :: mw

      INTEGER :: i
      INTEGER :: j
      INTEGER :: jj
      INTEGER :: k
      INTEGER :: l

      REAL :: delta
      REAL :: delta2
      REAL, DIMENSION(nlong,mnlat,nlev) :: global
      REAL, DIMENSION(mnlat,nmetlev)    :: intrp
      REAL, DIMENSION(mnlat,nlev)       :: global_2d
      REAL, DIMENSION(nlev)             :: eta_mid
      REAL, DIMENSION(mnlat)            :: lat_mid

      COMMON /SUM3_D/ global

! DEPENDS ON: sum3d
      CALL SUM3D(x)
      global_2d = SUM(global,DIM=1)*(mw/mair)/REAL(SIZE(global,DIM=1))

! Interpolate from mid point of STOCHEM grids,
!  so make appropriate ones
      eta_mid = (eta_stochem(1:nlev) + eta_stochem(0:nlev-1)) / 2.0
      lat_mid = (lat(1:mnlat) + lat(0:mnlat-1)) / 2.0

! Interpolate with height first
      DO k=1,nmetlev
        IF (eta_theta(k) < eta_mid(1)) THEN
          intrp(:,k) = global_2d(:,1)
        ELSE IF (eta_theta(k) > eta_mid(nlev)) THEN
          intrp(:,k) = global_2d(:,nlev)
        ELSE
          l = BISECT(eta_mid,eta_theta(k),nlev)
          delta = eta_mid(l+1) - eta_mid(l)
          delta2 = eta_theta(k) - eta_mid(l)
          intrp(:,k) = global_2d(:,l) + (global_2d(:,l+1) -             &
     &      global_2d(:,l)) * delta2 / delta
        END IF
      END DO

      DO k=1,nmetlev
        DO j=1,nlatpe
          jj = j + lobound - 1
          IF (latm_half(jj) < lat_mid(1)) THEN
            y(j,k) = intrp(1,k)
          ELSE IF (latm_half(jj) > lat_mid(mnlat)) THEN
            y(j,k) = intrp(mnlat,k)
          ELSE
            i = BISECT(lat_mid,latm_half(jj),mnlat)
            delta = dlat
            delta2 = latm_half(jj) - lat_mid(i)
            y(j,k) = intrp(i,k) + (intrp(i,k) - intrp(i+1,k)) *         &
     &        delta2 / delta
          END IF
        END DO
      END DO

      CONTAINS

      INTEGER FUNCTION BISECT(earray,z,ji)
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  5.5    14/02/03  Created. C.E. Johnson
!  6.1    06/08/04  Reformatted code. M.G. Sanderson
!
! Bisection method
!
      REAL, DIMENSION(*), INTENT(IN) :: earray
      REAL,    INTENT(IN) :: z
      INTEGER, INTENT(IN) :: ji

      INTEGER :: jm
      INTEGER :: jl
      INTEGER :: ju

      jl = 1
      ju = ji
      jm = (ju+jl)/2
      DO
        IF (ju-jl == 1) EXIT
        IF (z > earray(jm)) THEN
          jl = jm
        ELSE
          ju = jm
        END IF
        jm = (ju+jl)/2
      END DO

      BISECT = jm

      END FUNCTION BISECT

      END SUBROUTINE DATA2MET2D
#endif
