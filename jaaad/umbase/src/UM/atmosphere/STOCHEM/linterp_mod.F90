#if defined(A25_1A)
      MODULE LINTERP_MOD
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      PRIVATE

      INTERFACE LINTERP
         MODULE PROCEDURE LINTERP
      END INTERFACE ! LINTERP

      INTERFACE LINTERP_V
         MODULE PROCEDURE LINTERP_v
      END INTERFACE ! LINTERP_v

      PUBLIC LINTERP,LINTERP_v

      CONTAINS
      SUBROUTINE LINTERP_V(nfill,pos,x,val,lhalf,lvert,lu)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : DO 3D INTERPOLATIONS
!-
!-   Returned value  : LINTERP
!-   Inputs  : POS,X,LHALF,LVERT,LU
!-   Outputs :
!-   Controls:
!-
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  16/02/04   Created. K. Ketelsen
!  6.1  24/09/04   Reformatted for F90 compilation. M.G. Sanderson
!  6.2  31/01/06   Placed names of interface statement end blocks
!                  behind comments for portability. T. Edwards
!
!-
!VVV  V5.2  LINTERP 3/IX/01 - New dynamics version
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE HEIGHT_MOD
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTERFACE   ! interface module contains this function.
! DEPENDS ON: st_height
        INTEGER FUNCTION ST_HEIGHT(pos,eta_array)
          CHARACTER(*), INTENT(IN) :: eta_array
          REAL,         INTENT(IN) :: pos
        END FUNCTION ST_HEIGHT
! DEPENDS ON: getmetpoint
        REAL FUNCTION GETMETPOINT(pos,field,lhalf,lu)
          USE IN_STOCHEM_GRD
          REAL, DIMENSION(4),             INTENT(IN) :: pos
          REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: field
          LOGICAL,                        INTENT(IN) :: Lhalf  ! T for 1
          LOGICAL, INTENT(IN), OPTIONAL              :: LU     ! T for U
        END FUNCTION GETMETPOINT
      END INTERFACE

      INTEGER,              INTENT(IN) :: nfill
      REAL,DIMENSION(:),    INTENT(OUT):: val
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev),                         &
     &                      INTENT(IN) :: x        ! Interpolation array
      REAL, DIMENSION(4,nfill), INTENT(IN) :: pos      ! Position
      LOGICAL,              INTENT(IN) :: lhalf    ! Horiz=1/2 grids
      LOGICAL,              INTENT(IN) :: lvert    ! .T. for U grid
      LOGICAL, OPTIONAL,    INTENT(IN) :: lu       ! .T. for U grid

      REAL :: d3
      REAL :: deta
      REAL :: x5
      REAL :: x6
      INTEGER :: k,n1
      INTEGER :: i
      INTEGER,DIMENSION(nfill)   :: ka
      REAL,DIMENSION(nfill)      :: d3_a,deta_a

      IF (lvert) THEN
        CALL HEIGHT_ETA_RHO(pos,ka)
        DO i=1,nfill
          k = ka(i)
          IF (k==0) THEN
            deta=Eta_rho(k+1)
            d3=pos(3,i)
          ELSE
            deta=Eta_rho(k+1)-Eta_rho(k)
            d3=pos(3,i)-Eta_rho(k)
          END IF
          d3_a(i)   = d3
          deta_a(i) = deta
        END DO
      ELSE
        CALL HEIGHT_ETA_THETA(pos,ka)
        DO i=1,nfill
          k = ka(i)
          deta=Eta_theta(k+1)-Eta_theta(k)
          d3=pos(3,i)-Eta_theta(k)
          d3_a(i)   = d3
          deta_a(i) = deta
        END DO
      ENDIF

      IF (lhalf) THEN
        DO i=1,nfill
          k    = ka(i)
          d3   = d3_a(i)
          deta = deta_a(i)
! DEPENDS ON: getmetpoint
          x5 = GETMETPOINT(pos(:,i),x(:,:,k),lhalf)
! DEPENDS ON: getmetpoint
          x6 = GETMETPOINT(pos(:,i),x(:,:,k+1),lhalf)
          val(i)=x5+((x6-x5)/deta)*d3
        END DO
      ELSE
        DO i=1,nfill
          k    = ka(i)
          d3   = d3_a(i)
          deta = deta_a(i)
! DEPENDS ON: getmetpoint
          x5 = GETMETPOINT(pos(:,i),x(:,:,k),lhalf,lu)
! DEPENDS ON: getmetpoint
          x6 = GETMETPOINT(pos(:,i),x(:,:,k+1),lhalf,lu)
          val(i)=x5+((x6-x5)/deta)*d3
        END DO
      END IF

      END SUBROUTINE LINTERP_V

      REAL FUNCTION LINTERP(pos,x,lhalf,lvert,lu)
!-----------------------------------------------------------------------
!
!-   Purpose and Methods : DO 3D INTERPOLATIONS
!-
!-   Returned value  : LINTERP
!-   Inputs  : POS,X,LHALF,LVERT,LU
!-   Outputs :
!-   Controls:
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  3.0  09/12/93   Created. W.J. Collins
!  4.4  08/10/97   Converted to F90 and parallelised. W.J. Collins
!  5.5  28/11/01   Added a check for K==0 on rho g. Bill Collins
!  5.5  16/02/04   Refined as a module procedure and modified for
!                  vectorisation on SX6. K. Ketelsen
!  6.1  24/09/04   Reformatted for F90 compilation. M.G. Sanderson
!
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTERFACE   ! interface module contains this function.
! DEPENDS ON: st_height
        INTEGER FUNCTION ST_HEIGHT(pos,eta_array)
          CHARACTER(*), INTENT(IN) :: eta_array
          REAL,         INTENT(IN) :: pos
        END FUNCTION ST_HEIGHT
! DEPENDS ON: getmetpoint
        REAL FUNCTION GETMETPOINT(pos,field,lhalf,lu)
          USE IN_STOCHEM_GRD
          REAL, DIMENSION(4),             INTENT(IN) :: pos
          REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: field
          LOGICAL,                        INTENT(IN) :: Lhalf  ! T for 1
          LOGICAL, INTENT(IN), OPTIONAL              :: LU     ! T for U
        END FUNCTION GETMETPOINT
      END INTERFACE

      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev),                         &
     &                    INTENT(IN) :: x        ! Interpolation array
      REAL, DIMENSION(4), INTENT(IN) :: pos      ! Position
      LOGICAL,              INTENT(IN) :: lhalf    ! Horiz=1/2 grids
      LOGICAL,              INTENT(IN) :: lvert    ! .T. for U grid
      LOGICAL, OPTIONAL,  INTENT(IN) :: lu       ! .T. for U grid

      REAL :: d3
      REAL :: deta
      REAL :: x5
      REAL :: x6
      INTEGER :: k,n1

      IF (lvert) THEN
!kk     ST_HEIGHT expanded by hand
!kk        k=ST_HEIGHT(pos(3),'Eta_rho')
        n1=1
        IF (pos(3) > eta_rho(nmetlev)) THEN
          k=nmetlev
        ELSE
          DO
            IF (pos(3) <= eta_rho(n1)) EXIT
            n1=n1+1
          END DO
          k=n1-1
        END IF

        IF (k==0) THEN
          deta=Eta_rho(k+1)
          d3=pos(3)
        ELSE
          deta=Eta_rho(k+1)-Eta_rho(k)
          d3=pos(3)-Eta_rho(k)
        END IF
      ELSE
!kk        k=ST_HEIGHT(pos(3),'Eta_theta')
        n1=1
        DO
          IF (pos(3) <= eta_theta(n1)) EXIT
          n1=n1+1
        END DO
        k=n1-1

        deta=Eta_theta(k+1)-Eta_theta(k)
        d3=pos(3)-Eta_theta(k)
      END IF

      IF (lhalf) THEN
! DEPENDS ON: getmetpoint
        x5=GETMETPOINT(pos,x(:,:,k),lhalf)
! DEPENDS ON: getmetpoint
        x6=GETMETPOINT(pos,x(:,:,k+1),lhalf)
      ELSE
! DEPENDS ON: getmetpoint
        x5=GETMETPOINT(pos,x(:,:,k),lhalf,lu)
! DEPENDS ON: getmetpoint
        x6=GETMETPOINT(pos,x(:,:,k+1),lhalf,lu)
      END IF
      linterp=x5+((x6-x5)/deta)*d3

      END FUNCTION LINTERP

      END MODULE LINTERP_MOD
#endif
