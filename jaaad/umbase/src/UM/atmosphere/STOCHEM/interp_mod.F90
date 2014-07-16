#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      MODULE INTERP_MOD

      USE IN_STOCHEM_GRD
      USE HEIGHT_MOD

      PRIVATE

      INTERFACE INTERP_V
        MODULE PROCEDURE INTERP_V
        MODULE PROCEDURE INTERP_VU
      END INTERFACE ! INTERP_V

      PUBLIC INTERP_V

      CONTAINS

      SUBROUTINE INTERP_V(asize,todo,pos,x,val,lhalf,lvert)
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  16/02/04   Created. K. Ketelsen
!  6.1  07/09/04   Added '&' at end of lines for F90 compatability.
!                  M.G. Sanderson
!  6.2  31/01/06   Place names of interface blocks on end statements
!                  behind comments for portability. T. Edwards
!
      IMPLICIT  NONE
!
      INTEGER,                   INTENT(IN) :: asize
      LOGICAL, DIMENSION(asize), INTENT(IN) :: todo
      REAL, DIMENSION(asize),   INTENT(OUT) :: val
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev),                         &
     &                           INTENT(IN) :: x        ! Interpolation
      REAL, DIMENSION(4,asize),  INTENT(IN) :: pos      ! Position
      LOGICAL,            INTENT(IN) :: lhalf    ! Horiz=1/2 grids
      LOGICAL,            INTENT(IN) :: lvert    ! Vert=1/2 grids

      CALL INTERP_VU(asize,todo,pos,x,val,lhalf,lvert,.false.)

      RETURN

      END SUBROUTINE INTERP_V

      SUBROUTINE INTERP_VU(asize,todo,pos,x,val,lhalf,lvert,lu)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Do 3-D interpolation of wind and temperature
!-
!-   Returned value  :
!-   Inputs  : POS,X,LHALF,LVERT,LU
!            : U grid   LHALF=.F.,LVERT=.T.,LU=.T.
!            : V grid   LHALF=.F.,LVERT=.T.,LU=.F.
!            : W grid   LHALF=.T.,LVERT=.F.,LU=(not reqd.)
!-   Outputs :
!-   Controls:
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  16/02/04   Created. K. Ketelsen
!  6.1  07/09/04   Added '&' at end of lines for F90 compatability.
!                  M.G. Sanderson
!  6.2  20/04/05   Further improvements to avoid memory bank conflicts.
!                  R. Johanni.
!
!
!-
!VVV  V6.1  INTERP_VU 07/09/04
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------

      INTEGER,                   INTENT(IN)  :: asize
      LOGICAL, DIMENSION(asize), INTENT(IN)  :: todo
      REAL, DIMENSION(asize),    INTENT(OUT) :: val
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev),                         &
     &                           INTENT(IN) :: x   ! Interpolation array
      REAL, DIMENSION(4,asize),  INTENT(IN) :: pos ! Position
      LOGICAL,            INTENT(IN) :: lhalf    ! Horiz=1/2 grids
      LOGICAL,            INTENT(IN) :: lvert    ! Vert=1/2 grids
      LOGICAL,            INTENT(IN) :: lu       ! .T. for U grid, reqd.
                                                 ! only if LHALF=.F.
      INTEGER, PARAMETER :: u_grid=1
      INTEGER, PARAMETER :: v_grid=2
      INTEGER, PARAMETER :: w_grid=3
      INTEGER :: im   ! Longitude index on UM grid
      INTEGER :: jm   ! Latitude index on UM grid
      INTEGER :: grid
      INTEGER :: i
      INTEGER :: ii
      INTEGER :: klev
      INTEGER :: num
      INTEGER, DIMENSION(asize) :: ka
      INTEGER, DIMENSION(asize) :: idx

      REAL    :: d1   ! Longitudinal distance from grid point
      REAL    :: d2   ! Latitudinal distance from grid point
      REAL, DIMENSION(4) :: xv, yv
      REAL, DIMENSION(0:nmetlev,0:15) :: eta_x
      REAL, DIMENSION(asize) :: xx
!
! STOCHEM functions
      INTEGER :: IMET, JMET
      REAL    :: HINTERP, GETMETPOINT

      IF (lvert) THEN

! The horizontal wind components are on the rho grid which does not
! have a surface value, so we set an imaginary surface value of 0 for
! heights below the first eta_rho level.

        eta_x(0,0) = 0.0
        eta_x(1:nmetlev,0) = eta_rho

! HEIGHT_ETA_RHO is not converted to table lookup yet,
! therefore vectorised AHEIGHT is used.
        WHERE (todo)
          xx = pos(3,:)
        ELSEWHERE
          xx = 1.0
        END WHERE
! DEPENDS ON: aheight
        CALL AHEIGHT(asize,xx,'Eta_rho',ka)
      ELSE
        CALL HEIGHT_ETA_THETA(pos,todo,ka)
        eta_x(:,0) = eta_theta
      END IF

! The array eta_x contains only 21 elements but is referenced
! over and over in the loop below, leading to memory bank conflicts.
! As a little help for alleviating these conflicts we duplicate
! this array sometimes and the loop iterations are using one
! of the duplicates cyclically

      DO i = 1, 15
        eta_x(:,i) = eta_x(:,0)
      END DO

      IF (lhalf) THEN
        grid = w_grid
      ELSE IF (lu) THEN
        grid = u_grid
      ELSE
        grid = v_grid
      END IF

      DO i=1,asize
        IF (todo(i) .AND. ka(i) < nmetlev) THEN

! Function GETMETPOINT is inlined manually below,
! calling IMET and JMET only once
          IF (grid == u_grid) THEN          ! U
! DEPENDS ON: imet
            im = IMET(pos(1,i),.FALSE.)     ! 1
            d1 = pos(1,i) - longm(im+lnbound-1)
            IF (d1 < 0.0) d1 = d1 + 360.0
          ELSE                              ! V or W grid
! DEPENDS ON: imet
            im = IMET(pos(1,i),.TRUE.)      ! 1/2
            d1 = pos(1,i) - longm_half(im+lnbound-1)
            IF (d1 > 360.0) d1 = d1 - 360.0
          ENDIF
          IF (grid == v_grid) THEN          ! V grid
! DEPENDS ON: jmet
            jm = JMET(pos(2,i),.FALSE.)     ! 1
            d2 = pos(2,i) - latm(jm+lobound-1)
          ELSE                              ! U or W grid
! DEPENDS ON: jmet
            jm = JMET(pos(2,i),.TRUE.)      ! 1/2
            d2 = pos(2,i)- latm_half(jm+lobound-1)
          ENDIF

          ii = IAND(i,15)
          klev = MAX(ka(i)-1,0)             ! set to lowest available
! DEPENDS ON: hinterp
          yv(1) = HINTERP(x(:,:,klev),d1,d2,im,jm)
          xv(1) = eta_x(klev,ii)

          klev = ka(i)
! DEPENDS ON: hinterp
          yv(2) = HINTERP(x(:,:,klev),d1,d2,im,jm)
          xv(2) = eta_x(klev,ii)

          klev = ka(i)+1                    ! note that ka(i) < nmetlev
! DEPENDS ON: hinterp
          yv(3) = HINTERP(x(:,:,klev),d1,d2,im,jm)
          xv(3) = eta_x(klev,ii)

          klev = MIN(ka(i)+2,nmetlev)       ! set to highest available
! DEPENDS ON: hinterp
          yv(4) = HINTERP(x(:,:,klev),d1,d2,im,jm)
          xv(4) = eta_x(klev,ii)
! DEPENDS ON: cubfit
          CALL CUBFIT(xv,yv,pos(3,i),val(i))
        END IF
      END DO

! Care needed with cells with ka(i) == nmetlev:
! Assume met components do not vary above eta_rho top
! level, so no need for vertical interpolation.
! should only ever be true for eta_rho levels

      IF (lvert) THEN
! since only few of such cells exist, it is better to first
! gather these in an index array than having the IF in the loop
        num = 0
        DO i=1,asize
          IF (todo(i) .AND. ka(i) == nmetlev) THEN
            num = num + 1
            idx(num) = i
          END IF
        END DO
        DO i=1,num
          val(idx(i)) = &
! DEPENDS ON: getmetpoint
     &      GETMETPOINT(pos(:,idx(i)),x(:,:,nmetlev),lhalf,lu)
        END DO
      END IF

      END SUBROUTINE INTERP_VU

      END MODULE INTERP_MOD
#endif
