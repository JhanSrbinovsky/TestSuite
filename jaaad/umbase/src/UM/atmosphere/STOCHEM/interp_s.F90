#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      REAL FUNCTION INTERP_S(pos,x,lhalf,lvert,lu)
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
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.4    09/12/93  Created.  W.J. Collins
!  5.2    14/06/01  New Dynamics verion. C.E. Johnson
!  5.3    28/11/01  Split the LVERT and LHALF test for XR and EX1.
!                   Changed logic for J2==0 test.  W.J. Collins
!  5.5    22/01/04  Function ST_HEIGHT expanded by hand. K. Ketelsen
!  6.1    19/08/04  Changed name to INTERP_S to avoid clash with
!                   another UM variable. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson

!-
!VVV  V5.1  INTERP 10/IX/01 - Use GETMETPOINT etc
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTERFACE
! DEPENDS ON: getmetpoint
        REAL FUNCTION GETMETPOINT(pos,field,Lhalf,lu)
          USE IN_STOCHEM_GRD
          REAL, DIMENSION(4),             INTENT(IN) :: pos
          REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: field
          LOGICAL,                     INTENT(IN) :: Lhalf  ! T for 1/2
          LOGICAL,              OPTIONAL, INTENT(IN) :: lu
        END FUNCTION GETMETPOINT
! DEPENDS ON: st_height
        INTEGER FUNCTION ST_HEIGHT(pos,eta_array)
          REAL,         INTENT(IN) :: pos
          CHARACTER(*), INTENT(IN) :: eta_array
        END FUNCTION ST_HEIGHT
      END INTERFACE

      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev),                         &
     &                    INTENT(IN) :: x        ! Interpolation array
      REAL, DIMENSION(4), INTENT(IN) :: pos      ! Position
      LOGICAL,            INTENT(IN) :: lhalf    ! Horiz=1/2 grids
      LOGICAL,            INTENT(IN) :: lvert    ! Vert=1/2 grids
      LOGICAL, OPTIONAL,  INTENT(IN) :: lu       ! .T. for U grid, reqd.
!                                                ! only if LHALF=.F.

      INTEGER :: j
      INTEGER :: j2
      INTEGER :: k,n1

      REAL               :: x3
      REAL, DIMENSION(4) :: xr
      REAL, DIMENSION(4) :: ex1

      IF (lvert) THEN
!kk     Expand ST_HEIGHT by hand
!kk        k=ST_HEIGHT(pos(3),'Eta_rho')
        n1=1
        IF (POS(3) > Eta_rho(NMETLEV)) THEN
          k=NMETLEV
        ELSE
          DO
            IF(POS(3) <= Eta_rho(N1)) EXIT
            N1=N1+1
          END DO
          k=N1-1
        END IF
      ELSE
!kk        k = ST_HEIGHT(pos(3),'Eta_theta')
        n1=1
        DO
          IF (POS(3)<= eta_theta(n1)) EXIT
          N1=N1+1
        END DO
        k = n1-1
      ENDIF

! Assume met components do not vary above eta_rho top
! level, so no need for vertical interpolation.
      IF (k==nmetlev) THEN ! should only ever be true for eta_rho levels
        IF (.NOT. lhalf) THEN
! DEPENDS ON: getmetpoint
          INTERP_S = GETMETPOINT(pos,x(:,:,k),lhalf,lu)
        ELSE
! DEPENDS ON: getmetpoint
          INTERP_S = GETMETPOINT(pos,x(:,:,k),lhalf)
        END IF
      ELSE
! Do vertical cubic interpolation
! The horizontal wind components are on the rho grid which does not
! have a surface value, so we set an imaginary surface value of 0 for
! heights below the first Eta_rho level.
        DO j = 1,4
          j2 = k - 2 + j
          IF (j2 < 0) j2 = 0                  ! set to lowest available
          IF (j2 > nmetlev) j2 = nmetlev      ! set to highest available
          IF (.NOT. lhalf) THEN               ! V & U
! DEPENDS ON: getmetpoint
            xr(j) = GETMETPOINT(pos,x(:,:,j2),lhalf,lu)
          ELSE                                ! W and T
! DEPENDS ON: getmetpoint
            xr(j) = GETMETPOINT(pos,x(:,:,j2),lhalf)
          END IF
          IF (lvert) THEN
            IF (j2 == 0) THEN
              ex1(j) = 0.0
            ELSE
              ex1(j) = eta_rho(j2)
            END IF
          ELSE
            ex1(j) = eta_theta(j2)
          END IF
        END DO        ! J
! DEPENDS ON: cubfit
        CALL CUBFIT(ex1,xr,pos(3),x3)
        INTERP_S = x3
      END IF

      END FUNCTION INTERP_S
#endif
