#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine (Stochem) to convectively mix model cells.
!     Subroutine Interface:
!
      SUBROUTINE CLMIX2(pos,p,lnp,up,upd,orog,astep,seed3,nfill,        &
     &  z_top_of_model,first_constant_r_rho_level)
!
! Description:
!   This routine takes convective diagnostics supplied from the
!   atmosphere model and uses these to reposition model cells
!   in the vertical.  A randomised approach is used.
!
! Method:
!   Convective updraught mass flux and the updraught detrainment rate
!   are used to determine a range of destinations for upwards
!   movement. The cell enters into the updraught only if its random
!   number is greater than the convective fraction.  Its final
!   destination is also determined by random number.  All cells
!   in an active convective area subside by an amount given by
!   the change in updraught over adjacent levels.
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date                    Comment
!  5.5    28/02/03  Created.  W.J. Collins
!  5.5    23/02/04  Extensive changes to vectorise code. K. Ketelsen
!  6.1    20/08/04  STOCHEM upgrade.  C Johnson
!  6.2    06/04/05   Major changes for improved vectorisation.
!                                              R. Johanni
!
! Code description:
!   FORTRAN 90
!   This code is written to the programming standards of version 6
!   of UMDP3.
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      USE HEIGHT_MOD
      USE P2ETA_MOD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nfill
      INTEGER, INTENT(IN) :: first_constant_r_rho_level

      INTEGER, DIMENSION(ransize), INTENT(INOUT)            :: seed3
      REAL,                                      INTENT(IN) :: astep
      REAL,                                 INTENT(IN) :: z_top_of_model
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(INOUT) :: up
! We redefine UPD so it needs to be INOUT - a bit sloppy.
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(INOUT) :: upd
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev),  INTENT(IN) :: p
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev),  INTENT(IN) :: lnp
      REAL, DIMENSION(nlonpe,nlatpe),            INTENT(IN) :: orog
      REAL, DIMENSION(4,nclprc),              INTENT(INOUT) :: pos

      INTEGER :: im
      INTEGER :: jm
      INTEGER :: k
      INTEGER :: km
      INTEGER :: l
      INTEGER :: it
      INTEGER :: ktop

      REAL, PARAMETER :: cstep=450. ! 1/8 hour convective timestep
      REAL, PARAMETER :: minflux=1.0e-8 ! min non-zero massflux
      REAL :: deltap
      REAL :: influx
      REAL :: p1
      REAL :: p2
      REAL :: sub
      REAL :: up1
      REAL :: up2

      CHARACTER(LEN=72) :: cmessage

      INTEGER :: ntodo
      INTEGER :: nconv
      INTEGER, DIMENSION(nfill) :: im_a, jm_a, km_a, idxa
      INTEGER, DIMENSION(nfill) :: im_c, jm_c, km_c, idxc, contop, kk
      LOGICAL :: lerr
      LOGICAL, DIMENSION(nfill) :: lsub
      LOGICAL, DIMENSION(nlonpe,nlatpe) :: todo_flag
      REAL    :: zorog
      REAL, DIMENSION(nlonpe,nlatpe) :: eta2pnlev
      REAL, DIMENSION(nfill)    :: pos5
      REAL, DIMENSION(nfill)    :: confrac
      REAL, DIMENSION(nfill)    :: rnd_num
      REAL, DIMENSION(nfill)    :: aux
      REAL, DIMENSION(4,nfill)  :: xpos, kpos
      REAL, DIMENSION(nfill,nmetlev)   :: dest
      REAL, DIMENSION(nfill,0:nmetlev) :: sum_dest
      REAL :: gpos(4,nlonpe*nlatpe), pval(nlonpe*nlatpe)

! Precompute ETA2P value for NLEV
      ntodo = 0
      DO jm = 1, nlatpe-1   ! hinterp needs 2 values for interpolation
        DO im = 1, nlonpe-1 ! ditto
          ntodo = ntodo+1
          gpos(1,ntodo) = longm_half(im+lnbound-1)
          gpos(2,ntodo) = latm_half(jm+lobound-1)
          gpos(3,ntodo) = eta_stochem(nlev)
          gpos(4,ntodo) = 0.0
        END DO
      END DO
! DEPENDS ON: eta2p_v
      CALL ETA2P_V(ntodo,gpos,lnp,pval)
      ntodo = 0
      DO jm = 1, nlatpe-1
        DO im = 1, nlonpe-1
          ntodo = ntodo + 1
          eta2pnlev(im,jm) = pval(ntodo)
        END DO
      END DO

! Set up todo_flag ...
      todo_flag = .false.
      DO k=1,nmetlev
        DO jm = 1, nlatpe
          DO im = 1, nlonpe
            todo_flag(im,jm) =                                          &
     &        todo_flag(im,jm) .OR. (up(im,jm,k) > minflux)
          END DO
        END DO
      END DO

! ... and gather the values of all cells where there is something to do
      ntodo = 0
      DO l=1,nfill
! Met indicies - theta grid assumed.
! Use grid square values, not interpolated ones.
! DEPENDS ON: imet
        im = IMET(pos(1,l),.TRUE.)
! DEPENDS ON: jmet
        jm = JMET(pos(2,l),.TRUE.)
        IF (todo_flag(im,jm)) THEN
          ntodo = ntodo + 1
          xpos(1,ntodo) = pos(1,l)
          xpos(2,ntodo) = pos(2,l)
          xpos(3,ntodo) = pos(3,l)
          xpos(4,ntodo) = pos(4,l)
          kpos(1,ntodo) = longm_half(im+lnbound-1)
          kpos(2,ntodo) = latm_half(jm+lobound-1)
          im_a(ntodo) = im
          jm_a(ntodo) = jm
          idxa(ntodo) = l
        END IF
      END DO

! Top of models is now coincident, but may not always be
! DEPENDS ON: st_height
      ktop = ST_HEIGHT(eta_stochem(nlev),'Eta_theta')

!     WRITE(6,*) 'CLMIX Iterations: ',INT(astep/cstep)
      DO it = 1, INT(astep/cstep) ! several timesteps per advection step

! Check if there is a position out of bounds
        lerr = .false.
        DO l=1,ntodo
          IF (xpos(3,l) < eta_stochem(0) .OR.                           &
     &      xpos(3,l) > eta_stochem(nlev)) THEN
            lerr = .true.
          END IF
        END DO
        IF (lerr) THEN
          DO l=1,ntodo
            IF (xpos(3,l) < eta_stochem(0) .OR.                         &
     &        xpos(3,l) > eta_stochem(nlev)) THEN
              cmessage = 'Position out of bounds'
              WRITE(6,*) cmessage,'POS(3): ',xpos(3,l),' L: ',idxa(l)
! DEPENDS ON: ereport
              CALL EREPORT('CLMIX2',1,cmessage)
            END IF
          END DO
        END IF

! Check if there is any convection
        CALL HEIGHT_ETA_THETA(xpos,km_a(1:ntodo))
        kpos(3,1:ntodo) = xpos(3,1:ntodo)
! DEPENDS ON: eta2p_v
        CALL ETA2P_V(ntodo,kpos,lnp,pos5)
        DO l=1,ntodo
          im = im_a(l)
          jm = jm_a(l)
          km = km_a(l)
! get level thickness in Pa
          IF (km == ktop) THEN ! top level
            deltap = p(im,jm,km) - eta2pnlev(im,jm)
          ELSE
            deltap = p(im,jm,km) - p(im,jm,km+1)
          END IF
! in a convective plume ?
          IF (km < ktop .AND. up(im,jm,km+1) > 0.0) THEN
            IF (km == 0) THEN
              influx = up(im,jm,km+1)
            ELSE
              influx = up(im,jm,km+1) - up(im,jm,km) - upd(im,jm,km)
            END IF
          ELSE
            influx = 0.0
          END IF
! Fraction of column convecting.
          confrac(l) = influx * cstep / deltap
        END DO

! gather all cells which participate in convection
        lsub(1:ntodo) = .FALSE.
        nconv = 0
        DO l=1,ntodo
          IF(FRANV(seed3) < confrac(l)) THEN ! in updraught.
            nconv = nconv+1
            idxc(nconv) = l
            im_c(nconv) = im_a(l)
            jm_c(nconv) = jm_a(l)
            km_c(nconv) = km_a(l)
            lsub(l) = .TRUE.
          END IF
        END DO

! find top of the convection
        contop(1:nconv) = ktop
        DO k=1,ktop
          DO l=1,nconv
            im = im_c(l)
            jm = jm_c(l)
            km = km_c(l)
            IF (contop(l) == ktop .AND. k >= km+2 .AND. up(im,jm,k)==0) &
     &        contop(l) = k-1
          END DO
        END DO

! Detrain all mass at cloud top
!CDIR NODEP
        DO l=1,nconv
          im = im_c(l)
          jm = jm_c(l)
          k = contop(l)
          upd(im,jm,k) = -up(im,jm,k)
        END DO

! Normalised destination probabilites
        dest(1:nconv,:) = 0
        sum_dest(1:nconv,:) = 0
        DO k=1,ktop
          DO l=1,nconv
            im = im_c(l)
            jm = jm_c(l)
            km = km_c(l)
            IF (k >= km+1 .AND. k <= contop(l)) THEN
              dest(l,k) = (-upd(im,jm,k) / up(im,jm,k)) *               &
     &          (1.0-sum_dest(l,k-1))
              sum_dest(l,k) = sum_dest(l,k-1) + dest(l,k)
            END IF
          END DO
        END DO

! Randomly assign to a level according to probabilities
        DO l = 1, nconv
          rnd_num(l) = FRANV(seed3)
        ENDDO
        kk(1:nconv) = -1
        DO k=1,ktop
          DO l=1,nconv
! find random detrainment level
            IF(kk(l)<0 .AND. rnd_num(l) < sum_dest(l,k)) kk(l) = k
          END DO
        END DO
        IF (ANY(kk(1:nconv)<0)) THEN
          cmessage='**** LEVEL ERROR IN CLMIX2****'
          WRITE(6,*) cmessage
! DEPENDS ON: ereport
          CALL EREPORT('CLMIX2',1,cmessage)
        END IF

! Place cell randomly in level in pressure
        DO l=1,nconv
          im = im_c(l)
          jm = jm_c(l)
          IF (kk(l) == ktop) THEN
            p2 = eta2pnlev(im,jm)
          ELSE
            p2 = p(im,jm,kk(l)+1)
          END IF
          p1 = p(im,jm,kk(l))
          pos5(idxc(l)) =                                               &
     &     p1+(p2-p1)*(rnd_num(l)-sum_dest(l,kk(l)-1)) / dest(l,kk(l))
        END DO

        CALL P2ETA_V(ntodo,lsub,pos5,kpos,p,xpos)

! All  cells subside by a fraction of a level given by UP
        CALL HEIGHT_ETA_THETA(xpos,km_a(1:ntodo))
        DO l=1,ntodo
          im = im_a(l)
          jm = jm_a(l)
          km = km_a(l)
          IF (km == 0) THEN ! bottom level
            p1 = p(im,jm,0)
            up1 = 0.0
          ELSE
            p1 = p(im,jm,km)
            up1 = up(im,jm,km)
          END IF
          IF (km == ktop) THEN ! top level
            p2 = eta2pnlev(im,jm)
            up2 = 0.0
          ELSE
            p2 = p(im,jm,km+1)
            up2 = up(im,jm,km+1)
          END IF

! Linearly interpolate subsidence fluxes
          sub = (up1*(pos5(l)-p2) + up2*(p1-pos5(l))) / (p1-p2)

! Do subsidence in pressure rather than eta for accuracy.
! Check new pressure is not greater than p0.
          lsub(l) = .FALSE.
          IF (sub > minflux) THEN
            lsub(l) = .TRUE.
            aux(l) = pos5(l) + sub * cstep
          END IF
        END DO

        CALL P2ETA_V(ntodo,lsub,aux,kpos,p,xpos)

      END DO

      DO l=1,ntodo
! DEPENDS ON: getmetpoint
        zorog = GETMETPOINT(xpos(:,l),orog,.TRUE.)
! DEPENDS ON: etator
        xpos(4,l) = ETATOR(xpos(3,l),zorog,z_top_of_model,              &
     &    first_constant_r_rho_level)
      END DO

      DO l=1,ntodo
        pos(3,idxa(l)) = xpos(3,l)
        pos(4,idxa(l)) = xpos(4,l)
      END DO

      END SUBROUTINE CLMIX2
#endif
