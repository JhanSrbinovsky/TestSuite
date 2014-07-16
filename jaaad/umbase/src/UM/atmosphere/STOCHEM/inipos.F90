#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE INIPOS(pos,cellno,nfill,orog,z_top_of_model,           &
     &  first_constant_r_rho_level)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : SET UP INITIAL CELL POSITIONS
!-
!-   Inputs  : PROCMAP, LAT (from module)
!-   Outputs : POS,NFILL,CELLNO
!-   Controls:
!-
!
! History:
! Version   Date                    Comment
!  3.4    09/12/93  Created.  W.J. Collins
!  5.0    09/07/01  New Dynamics version. C.E. Johnson
!  5.5    06/01/04  Vectorised code. M.G. Sanderson
!  6.1    20/10/04  Now uses system initial random number seed.
!                   M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson

!VVV  V5.0  INIPOS 9/VII/01 - ND version
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: first_constant_r_rho_level
      REAL,                           INTENT(IN) :: z_top_of_model
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: orog
      REAL, DIMENSION(4,nclprc),     INTENT(OUT) :: pos
      INTEGER, DIMENSION(nclprc),    INTENT(OUT) :: cellno
      INTEGER,                       INTENT(OUT) :: nfill

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: k0
      INTEGER :: n
      INTEGER :: im
      INTEGER :: jm
      INTEGER :: mcell
      INTEGER, DIMENSION(0:mnlat-1,nlev) :: nn ! 2D cell distribution

      REAL :: xt
      REAL :: xgrid
      REAL :: zorog
      REAL :: deta
      REAL :: xran
      REAL :: t
      REAL :: reus    ! Temporary store
      REAL :: dz
      REAL, DIMENSION(mnlat,nlev) :: mfract
      REAL, DIMENSION(0:nlev)     :: press
      REAL, DIMENSION(0:nlev)     :: tstd =                             &
                                             ! standard atmosphere
     &  (/287.5,284.6,281.6,278.7,275.5,272.4,268.5,265.0,265.0,260.9,  &
     &    256.5,252.0,247.0,241.7,235.5,228.7,220.0,216.6,216.6,216.6,  &
     &    217.2/)
      COMPLEX :: seed  ! Random number seed

! Set up initial random number seed. Values from call to FRANDI()
      seed = (2.224615814526098E+41,-2.001516614211545E+118)

! Calculate pressure on stochem levels
      press(0) = pstar
      DO k = 1, nlev
        t = (tstd(k) + tstd(k-1)) / 2.0
        dz = (eta_stochem(k) - eta_stochem(k-1)) * z_top_of_model
        press(k) = press(k-1) * EXP(-1.0*g*dz/(rmol*t/mair))
      END DO

! Calculate no. of cells in each band
      DO j = 1, mnlat
        mfract(j,nlev) =                                                &
     &   SIN((lat(j)-90.0)*pi_over_180) -                               &
     &   SIN((lat(j-1)-90.0)*pi_over_180)
      END DO
      reus = 1.0 / (2.0 * (press(0) - press(nlev)))
      DO k = 1, nlev
        mfract(:,k) = mfract(:,nlev) * (press(k-1)-press(k)) * reus
      END DO
      nn = INT(ncell*mfract)
      n = SUM(nn)
      xt = SUM(mfract)
      WRITE(6,*) ' *** INIPOS: Initially, N = ',n,sum(nn)

! Add in extra cells into random positions:
      k0 = n
      IF (k0 < ncell) THEN
        DO i = k0+1, ncell
          xran = FRAND(seed)
          j = INT(xran*mnlat)
          xran = FRAND(seed)
          k = 1 + INT(xran*nlev)
          nn(j,k) = nn(j,k) + 1
          n = n + 1
        END DO
        WRITE(6,*) ' *** INIPOS: Added cells, N = ',n,sum(nn)

! Subtract cells if there are too many:
      ELSE IF (k0 > ncell) THEN
        DO i = ncell+1, k0
          xran = FRAND(seed)
          j = INT(xran*mnlat)
          xran = FRAND(seed)
          k = 1 + INT(xran*nlev)
          DO
            IF (nn(j,k) > 0) EXIT
            j = j + 1
            IF (j > mnlat-1) THEN
              j = 0
              k = k + 1
            END IF
          END DO
          nn(j,k) = nn(j,k) - 1
          n = n - 1
        END DO
        WRITE(6,*) ' *** INIPOS: Removed cells, N = ',n,sum(nn)
      END IF

      mcell = n
      IF (mcell /= ncell) WRITE(6,*) '**** MCELL, NCELL=', mcell, ncell

! Write out initial conditions and allocation tables
      WRITE(6,'(A6,I3,/,A17,I8,A17,F8.2,/,A30,/)')                      &
     &  'NLEV: ',nlev,                                                  &
     &  'TOTAL NO. CELLS: ',n, ' TOTAL FRACTION: ', xt,                 &
     &  'ALLOCATION OF CELLS IN 1-D:'
      DO k = nlev, 1, -1
! Next two lines write out 1D (height) distribution of air parcels
        WRITE(6,'(A7,I6,A17,I6)')                                       &
     &    'LEVEL: ',k,' TOTAL IN LEVEL: ',SUM(nn(:,k))
! Next two lines write out 2D (lat-height) distribution
!       WRITE(6,'(A7,I6,A17,I6,/,6(12I6,/))')
!     &   'LEVEL: ',k,' TOTAL IN LEVEL: ',SUM(nn(:,k)),nn(0:mnlat-1,k)
      END DO

! Place cells in Eulerian grid squares with random component
      n = 1
      pos = -999.0
      DO j = 0, mnlat-1
        DO k = 1, nlev
          deta = eta_stochem(k) - eta_stochem(k-1)
          IF (nn(j,k) > 0)  THEN
            xgrid = 360.0 / nn(j,k)
            DO i = 1, nn(j,k)
              xran = FRAND(seed)
              pos(1,n) = xgrid*REAL(i-1) + xran*xgrid
              xran = FRAND(seed)
              pos(2,n) = lat(j) + xran*dlat
              xran = FRAND(seed)
              pos(3,n) = eta_stochem(k-1) + xran*deta

! Set up CELLNO
              cellno(n) = SUM(nn(0:j-1,:)) + SUM(nn(j,1:k-1)) + i
              im = 1 + INT(pos(1,n)/dlong)
              jm = 1 + INT(pos(2,n)/dlat)

! If cell is on this PE set R, then count it in.
              IF (procmap(im,jm) == mype) THEN
! Set R=POS(4,N)
! DEPENDS ON: getmetpoint
                zorog = GETMETPOINT(pos(:,n),orog,.true.)  ! 1/2 grid
! DEPENDS ON: etator
                pos(4,n) = ETATOR(pos(3,n),zorog,z_top_of_model,        &
     &            first_constant_r_rho_level)
                n = n + 1
              END IF
            END DO
          END IF
        END DO
      END DO

      nfill = n - 1   ! Number of air parcels on this PE

      END SUBROUTINE INIPOS
#endif
