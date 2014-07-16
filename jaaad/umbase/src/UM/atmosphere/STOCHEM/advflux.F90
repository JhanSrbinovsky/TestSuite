#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE ADVFLUX(pos,tropz,xx,flist,totflu,cellflux,nflux,      &
     & posold,nfill,lnp)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Adds up 3D fluxes
!-
!-   Inputs  : POS,POSOLD,FLIST,NFLUX,NFILL,LNP,TROPZ
!-   Outputs : TOTFLU,CELLFLUX
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.5    18/11/98  Created.  W.J. Collins
!  5.2    12/09/01  More accurate tropopause calculation. C.E. Johnson
!  6.1    21/07/04  Improved vectorisation. M.G. Sanderson
!
!  6.2    06/04/05   Further vectorisation. R. Johanni
!-
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_OUT
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER,                       INTENT(IN) :: nflux
      INTEGER,                       INTENT(IN) :: nfill
      INTEGER, DIMENSION(2,numflux), intent(in) :: flist

      REAL, DIMENSION(4,nclprc),                INTENT(IN) :: pos
      REAL, DIMENSION(4,nclprc),                INTENT(IN) :: posold
      REAL, DIMENSION(nc,nclprc),               INTENT(IN) :: xx
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN) :: tropz
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: lnp
      REAL, DIMENSION(num3dflux_dim,nclprc),   INTENT(OUT) :: cellflux
      REAL, DIMENSION(numflux),              INTENT(INOUT) :: totflu

      INTEGER :: j
      INTEGER :: jj
      INTEGER :: n
      INTEGER :: ntodo
      INTEGER, DIMENSION(nfill) :: ipos3
      INTEGER, DIMENSION(nfill) :: ipos3old
      INTEGER, DIMENSION(nfill) :: idx
      INTEGER, DIMENSION(nfill) :: flag
      REAL, DIMENSION(nfill) :: val
      REAL, DIMENSION(nfill) :: valold

      DO j = 1, num3dflux_dim
        DO n = 1, nclprc
          cellflux(j,n) = 0.0
        END DO
      END DO

! Section to diagnose advection fluxes.
      DO j=1,nfill
        flag(j) = 0
! If now above tropopause
        IF (pos(4,j) - earth_radius >                                   &
! DEPENDS ON: getmetpoint
     &    GETMETPOINT(pos(:,j),tropz,.TRUE.)) THEN
! and were not before
          IF (posold(4,j)-earth_radius <                                &
! DEPENDS ON: getmetpoint
     &       GETMETPOINT(posold(:,j),tropz,.TRUE.)) THEN
              flag(j) = -1
          END IF
        ELSE
! If now not in stratosphere and were before
          IF (posold(4,j)-earth_radius >                                &
! DEPENDS ON: getmetpoint
     &      GETMETPOINT(posold(:,j),tropz,.TRUE.)) THEN
            flag(j) = 1
          END IF
        END IF
      END DO

! Now sort out the cells where fluxes are ...

      ntodo = 0
      DO j=1,nfill
        IF (flag(j) /= 0) THEN
          ntodo = ntodo + 1
          idx(ntodo) = j
        END IF
      END DO

! ... and enter these into totflu and cellflux

      DO jj=1,ntodo
        j = idx(jj)
        IF (flag(j) < 0) THEN
          DO n=1,nflux
            IF (INT(flist(1,n)/100) == 7) THEN    ! STE fluxes
              totflu(n) = totflu(n) - xx(flist(1,n)-700,j)
              IF (flist(2,n) > 0) THEN            ! 3D Fluxes
                cellflux(flist(2,n),j) = -xx(flist(1,n)-700,j)
              END IF
            END IF
          END DO
        ELSE
          DO n=1,nflux
            IF (INT(flist(1,n)/100) == 7) THEN    ! STE fluxes
              totflu(n) = totflu(n) + xx(flist(1,n)-700,j)
              IF (flist(2,n) > 0) THEN            ! 3D Fluxes
                cellflux(flist(2,n),j) = xx(flist(1,n)-700,j)
              END IF
            END IF
          END DO
        END IF
      END DO

! Now calculate fluxes through pressure levels
! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,pos,lnp,val)
! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,posold,lnp,valold)

! Go through all cells and flag those where there is further
! work to be done
      DO j=1,nfill
        ipos3(j)    = INT(val(j)/10000.)     ! 10000 Pa spacing
        ipos3old(j) = INT(valold(j)/10000.)
        flag(j) = 0
        IF (ipos3(j) < 10 .AND. ipos3old(j) < 10) THEN
          IF (ipos3(j) < ipos3old(j)) THEN ! Upward fluxes (-ve)
            flag(j) = -1
          ELSE IF (ipos3(j) > ipos3old(j)) THEN ! Downward fluxes (+ve)
            flag(j) = 1
          END IF
        END IF
      END DO

! Now sort out the cells where fluxes are ...

      ntodo = 0
      DO j=1,nfill
        IF (flag(j) /= 0) THEN
          ntodo = ntodo + 1
          idx(ntodo) = j
        END IF
      END DO

! ... and enter these into totflu and cellflux

      DO jj=1,ntodo
        j = idx(jj)
        IF (flag(j) < 0) THEN
          DO n=1,nflux
            IF (INT(flist(1,n) / 100) == 7+ipos3old(j)) THEN ! Advection
              totflu(n) = totflu(n) -                                   &
     &          xx(flist(1,n)-700-100*ipos3old(j),j)
              IF (flist(2,n) > 0) THEN                 ! 3D Fluxes
                cellflux(flist(2,n),j) =                                &
     &            -xx(flist(1,n)-700-100*ipos3old(j),j)
              END IF
            END IF
          END DO
        ELSE
          DO n=1,nflux
            IF (INT(flist(1,n)/100) == 7+ipos3(j)) THEN ! Advection flux
              totflu(n) = totflu(n) +                                   &
     &           xx(flist(1,n)-700-100*ipos3(j),j)
              IF (flist(2,n) > 0) THEN                 ! 3D Fluxes
                cellflux(flist(2,n),j) =                                &
     &            xx(flist(1,n)-700-100*ipos3(j),j)
              END IF
            END IF
          END DO
        END IF
      END DO

      END SUBROUTINE ADVFLUX
#endif
