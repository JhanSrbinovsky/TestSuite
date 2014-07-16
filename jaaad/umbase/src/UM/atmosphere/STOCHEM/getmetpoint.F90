#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      REAL FUNCTION GETMETPOINT(pos,field,lhalf,lu)
! ----------------------------------------------------------------------
! Purpose:
! To interpolate a met field to a point
! Restricted to data on half grid points at present

! Method:
! Bi-linear

! Original Programmer: Colin Johnson
!
! Current code owner: Colin Johnson

! History:
! Date        Version     Comment
! -------     -------     -------
!  7/9/01     1.0         Original              Colin Johnson

!VVV  V1.1  GETMETPOINT  7/9/01 - Original
! ----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
! ----------------------------------------------------------------------
      INTERFACE
! DEPENDS ON: hinterp
        REAL FUNCTION HINTERP(x,d1,d2,i,l)
          USE IN_STOCHEM_GRD
          REAL, DIMENSION(lnbound:lnbound+nlonpe-1,                     &
     &                    lobound:lobound+nlatpe-1),                    &
     &                 INTENT(IN) :: x        ! Interpolation array
          REAL,        INTENT(IN) :: d1       ! distance from long i
          REAL,        INTENT(IN) :: d2       ! distance from lat l
          INTEGER,     INTENT(IN) :: i        ! longitude index
          INTEGER,     INTENT(IN) :: l        ! latitude index
        END FUNCTION HINTERP

! DEPENDS ON: imet
        INTEGER FUNCTION IMET(pos1,lhalf)
          USE IN_STOCHEM_GRD
          REAL,    INTENT(IN) :: pos1       ! longitude
          LOGICAL, INTENT(IN) :: lhalf      ! half grid or not
        END FUNCTION IMET

! DEPENDS ON: jmet
        INTEGER FUNCTION JMET(pos2,lhalf)
          USE IN_STOCHEM_GRD
          REAL,    INTENT(IN) :: pos2       ! latitude
          LOGICAL, INTENT(IN) :: lhalf      ! half grid or not
        END FUNCTION JMET
      END INTERFACE

      REAL, DIMENSION(4),             INTENT(IN) :: pos
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: field
      LOGICAL,                        INTENT(IN) :: Lhalf  ! T for 1/2
      LOGICAL,              OPTIONAL, INTENT(IN) :: lu     ! T for U

      INTEGER            :: im   ! Longitude index on UM grid
      INTEGER            :: jm   ! Latitude index on UM grid
      INTEGER, PARAMETER :: u_grid=1
      INTEGER, PARAMETER :: v_grid=2
      INTEGER, PARAMETER :: w_grid=3
      INTEGER            :: grid


      REAL    :: d1   ! Longitudinal distance from grid point
      REAL    :: d2   ! Latitudinal distance from grid point


      IF(lhalf) THEN
        grid=w_grid
      ELSEIF(lu) THEN
        grid=u_grid
      ELSE
        grid=v_grid
      ENDIF

      IF(grid==u_grid) THEN             ! U
! DEPENDS ON: imet
        im=IMET(pos(1),.FALSE.)   ! 1
        d1=pos(1)-longm(im+lnbound-1)
        IF(d1<0.0) d1=d1+360.0
      ELSE                              ! V or W grid
! DEPENDS ON: imet
        im=IMET(pos(1),.TRUE.)    ! 1/2
        d1=pos(1)-longm_half(im+lnbound-1)
        IF(d1>360.) d1=d1-360.0
      ENDIF
      IF (grid==v_grid) THEN            ! V grid
! DEPENDS ON: jmet
        jm=JMET(pos(2),.FALSE.)   ! 1
        d2=pos(2)-latm(jm+lobound-1)
      ELSE                              ! U or W grid
! DEPENDS ON: jmet
        jm=JMET(pos(2),.TRUE.)    ! 1/2
        d2=pos(2)- latm_half(jm+lobound-1)
      ENDIF
! Polar stuff now taken into account in halos.
! DEPENDS ON: hinterp
      GETMETPOINT=HINTERP(field,d1,d2,im,jm)

      END FUNCTION GETMETPOINT
#endif
