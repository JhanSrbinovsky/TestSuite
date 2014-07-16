#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE STRATCALC3(xx,m,em,pos,ipos,nfill,tropz,orog,o3um,     &
     &  astep)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : CALCULATE O3 & HNO3 influx from stratosphere
!-
!-   Inputs  : POS,TROPZ,OROG,O3_LS
!-
!-   Outputs : Set TROP & STRAT tracers in XX,
!-             Set EM for O3, STRAT_O3,HNO3, NO2, STRAT & TROP.
!
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.2    22/02/00  Created. W.J. Collins
!  5.5    16/07/02  Uses TROPZ instead of TROP
!  6.1    21/10/04  Refromatted code. M.G. Sanderson
!-
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER,                          INTENT(IN) :: nfill
      INTEGER, DIMENSION(3,nclprc),     INTENT(IN) :: ipos

      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: orog
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: tropz
      REAL, DIMENSION(4,nclprc),        INTENT(IN) :: pos
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(IN) :: o3um
      REAL,                                INTENT(IN) :: astep
      REAL, DIMENSION(nclprc),             INTENT(IN) :: m
      REAL, DIMENSION(nc,nclprc),       INTENT(INOUT) :: xx
      REAL, DIMENSION(nc+2,nclprc),     INTENT(INOUT) :: em

      INTEGER :: i,j,k,l
      REAL :: o3strat
      REAL :: noystrat
      REAL :: deltao3
      REAL :: deltanoy
      REAL :: deltano2
      REAL :: o3top
      REAL :: o3bot
! Ratio of kg[N]/kg[O3] from Murphy and Fahey 1994
      REAL, PARAMETER :: noy_frac=1.0/1000.0
! 20 day relaxation time
      REAL, PARAMETER :: tau20=20.0*daysec
! 2 day relaxation time
      REAL, PARAMETER :: tau2=2.0*daysec
      REAL    :: z1
      REAL    :: z2
      REAL    :: z3
      REAL    :: dz
      REAL    :: ztop
      REAL    :: zbot
      REAL    :: zorog
      REAL, DIMENSION(4) :: pos2
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev) :: o3dash

      em(i_o3,:) = 0.0
      em(i_hno3,:) = 0.0
      em(i_no2,:) = 0.0
      em(i_trop,:) = 0.0
      em(i_strat,:) = 0.0
      o3dash(:,:,0) = o3um(:,:,1)
      o3dash(:,:,1:nmetlev) = o3um

      DO j = 1,nfill
        pos2(:) = pos(:,j)
! DEPENDS ON: getmetpoint
        zorog = GETMETPOINT(pos2,orog,.TRUE.)
        IF (pos(4,j)-zorog-earth_radius >                               &
! DEPENDS ON: getmetpoint
     &    GETMETPOINT(pos2,tropz,.TRUE.)) THEN

! Use ozone ancillary data above tropopause:

! DEPENDS ON: interp_s
          o3strat = INTERP_S(pos(:,j),o3dash,.TRUE.,.FALSE.)

          noystrat = o3strat*noy_frac*(mo3/mnit)
          deltao3 = o3strat - xx(i_o3,j)
          deltanoy = noystrat - xx(i_hno3,j)
!          deltano2 = noystrat - xx(i_no2,j)
          em(i_o3,j) = deltao3/tau20
          em(i_hno3,j) = deltanoy/tau20
!          em(i_no2,j) = deltano2/tau20
          xx(i_strat,j) = 1.0
          em(i_trop,j) = -xx(i_trop,j)/tau2
        ELSE
! Below troposphere
          xx(i_trop,j) = 1.0
          em(i_strat,j) = -xx(i_strat,j)/tau2
        END IF
      END DO

      END SUBROUTINE STRATCALC3
#endif
