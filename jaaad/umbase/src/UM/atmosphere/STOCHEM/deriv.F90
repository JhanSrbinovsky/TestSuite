#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE DERIV(nflux,nfill,flist,astep,time,tl,q,rhl,md,liq,dj, &
     &  dryd,wetd,emit,xx,totflu,cellflux)
!-----------------------------------------------------------------------
!
! Current Owner of Code: M.G. Sanderson
!
! History:
! Version   Date                    Comment
!  3.4    01/09/94  Created.  C.E. Johnson
!  4.4    11/10/96  O+NO2 added. N2O5+H2O deleted but NAER added.
!                   HO2+NO3 added. HO2NO2 & CH3OH added, acetone,
!                   peroxides and DMS added. R.G. Derwent.
!  5.0    11/11/98  New species O3(S), CO(E). H2O in Y() array.
!                   W.J. Collins
!  5.1    04/01/00  Ozone budget terms added to flux array. C.E. Johnson
!  5.5    10/04/04  Code vectorised. K. Ketelsen
!  5.5    12/05/04  Reactions RO2IP1/RO2IP2 + CH3O2 removed (R242,R243).
!                   M.G. Sanderson
!  6.1    21/07/04  Reactions stopped in stratosphere. M.G. Sanderson
!  6.1    28/02/06  Detects parcels in clouds for faster execution.
!                   M.G. Sanderson
!  6.2    02/03/06  Loops over parcels in chunks instead of entire
!                   arrays.  M.G. Sanderson
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_OUT
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nflux                      ! No. 3D fluxes
      INTEGER, INTENT(IN) :: nfill                      ! No. cells on P
      INTEGER, DIMENSION(2,numflux), INTENT(IN) :: flist
      REAL, INTENT(IN) :: astep                         ! Timestep
      REAL, INTENT(IN) :: time                          ! Current time
      REAL, DIMENSION(nclprc), INTENT(IN) :: tl         ! Temperature of
      REAL, DIMENSION(nclprc), INTENT(IN) :: q          ! Water vapour c
      REAL, DIMENSION(nclprc), INTENT(IN) :: rhl        ! Relative humid
      REAL, DIMENSION(nclprc), INTENT(IN) :: md         ! Molecular dens
      REAL, DIMENSION(nclprc), INTENT(IN) :: liq        ! Liquid water c
      REAL, DIMENSION(ndj,nclprc), INTENT(IN) :: dj     ! Photolysis rat
! dd, dw, em shouldn't really be INOUT as we ought to avoid redefining t
      REAL, DIMENSION(nc,nclprc), INTENT(INOUT) :: dryd ! Dry deposition
      REAL, DIMENSION(nc,nclprc), INTENT(INOUT) :: wetd ! Wet deposition
      REAL, DIMENSION(nc+2,nclprc), INTENT(INOUT) :: emit ! Emission rat
      REAL, DIMENSION(nc,nclprc), INTENT(INOUT) :: xx   ! Species concen
      REAL, DIMENSION(numflux), INTENT(INOUT) :: totflu ! Global total f
      REAL, DIMENSION(num3dflux_dim,nclprc), INTENT(INOUT) :: cellflux

      INTEGER, PARAMETER :: nit=8  ! Number of iterations for Backward E
      INTEGER :: i                 ! Loop count over iterations
      INTEGER :: j                 ! Loop count over air parcels
      INTEGER :: k                 ! Loop count
      INTEGER :: kk                ! Loop count
      INTEGER :: j0                ! Lower index of air parcels (steps o
      INTEGER :: j1                ! Upper index of air parcel
      INTEGER :: asize             ! Max index for local arrays
      INTEGER :: k1                ! Set to flist(1,..)
      INTEGER :: k2                ! Set to flist(2,..)
      INTEGER :: ngtlcrit          ! Number of air parcels where ll > l_
      INTEGER, DIMENSION(chunk) :: indx

      REAL, PARAMETER :: l_crit = 1.0e-10   ! l water / l air
      REAL, PARAMETER :: cvf = 0.5          ! Cloud volume fraction
      REAL, PARAMETER :: one_over_cvf = 1.0 / cvf ! 1/(Cloud volume frac
      REAL :: small_conc = 4.0e-18
      REAL :: one_ppbv = 1.0e-9             ! 1 ppbv as a vmr
      REAL :: liq_const    ! Factor needed to convert liq to vol water/v
      REAL :: molec2mol    ! Converts molecules cm-3 to mol/l
      REAL :: mol2molec    ! Converts mol/l to molecules cm-3
      REAL :: molec2atm    ! Converts molecules cm-3 to pressure in atm
      REAL :: dts          ! Chemistry time step (s)
      REAL :: ctime        ! Time at beginning of step
      REAL :: tmax         ! Time at end of step
      REAL :: dtsm         ! dts / m
      REAL :: ee           ! Fraction of SO2 oxidised per timestep
      REAL :: pt, lt       ! Prod/Loss of OH, O(1D) and O(3P)
      real :: mtemp
      REAL, DIMENSION(chunk) :: f
      REAL, DIMENSION(chunk,nr) :: rc  ! Rate constants, cm, molecule, s
      REAL, DIMENSION(chunk,ndj) :: kj ! Photolysis rates for cells (s-1
      REAL, DIMENSION(chunk,nc) :: y   ! Current concentrations
      REAL, DIMENSION(chunk,nc) :: yp  ! Concentrations at previous step
      REAL, DIMENSION(chunk,6) :: kh   ! Henry's law constants (mol / l.
      REAL, DIMENSION(chunk,6) :: ke   ! Equilibrium constants
      REAL, DIMENSION(chunk,10) :: laq
      REAL, DIMENSION(chunk,5) :: faq
      REAL, DIMENSION(chunk,12) :: cs
      REAL, DIMENSION(chunk) :: tc
      REAL, DIMENSION(chunk) :: ql
      REAL, DIMENSION(chunk) :: m
      REAL, DIMENSION(chunk,nc) :: dd
      REAL, DIMENSION(chunk,nc) :: dw
      REAL, DIMENSION(chunk,nc+2) :: em
      REAL, DIMENSION(chunk,700) :: flux
      REAL, DIMENSION(chunk) :: p      ! Production rate of a species
      REAL, DIMENSION(chunk) :: l      ! Loss rate of a species
      REAL, DIMENSION(chunk) :: l1,l2,l3,r1,r2,ll
      REAL, DIMENSION(chunk) :: cth2o2,cthno3,ctnh3,ctnh4,so4,cpso4,    &
     &  cph2o2,ctso2,cpso2,cpnh3,cpnh4,c1h2o2,c1so2,c1so4,c1nh3,        &
     &  c1nh4,reus,fso2,hp,ph
      REAL, DIMENSION(chunk) :: dtsm_a
      LOGICAL, DIMENSION(chunk) :: in_strat

! some of the local variables are used before being defined,
! so set them all to 0 to avoid strange side effects
! (unless explicitly set to 0 somewhere else)

      f=0; rc=0; y=0; yp=0; kh=0; ke=0; laq=0; faq=0; p=0; l=0; kj=0
      l1=0; l2=0; l3=0; r1=0; r2=0
      cth2o2=0; cthno3=0; ctnh3=0; ctnh4=0; so4=0; cpso4=0; cph2o2=0
      ctso2=0; cpso2=0; cpnh3=0; cpnh4=0; c1h2o2=0; c1so2=0; c1so4=0
      c1nh3=0; c1nh4=0; reus=0; fso2=0; hp=0; ph=0; dtsm_a=0

! Set up constants
      liq_const = mair / (na * rho_h2o)   ! Factor to convert liq to vol
      molec2mol = 1.0e3 / na              ! molecules/cm3 to mol/l(air).
      mol2molec = 1.0 / molec2mol         ! mol/l(air) to molecules/cm3.
      molec2atm = rgc * 1.0e6 / (na * pstar) ! molecules/cm3 to atm

! Check for large ozone depletion by dry deposition.
! Do we still need this check ???
      DO j0 = 1, nfill
        IF (dryd(i_o3,j0) * astep > 0.2) dryd(i_o3,j0) = 0.2 / astep
      END DO

! Avoid ozone deposition if ozone < 1 ppb.
      DO j0 = 1, nfill
        IF (xx(i_o3,j0) < one_ppbv) dryd(i_o3,j0) = 0.0
      END DO

! Start of main loop
      DO j0 = 1, nfill, chunk  ! Loop over all air parcels with a step o
                               ! chunk
        j1 = MIN(nfill,j0+chunkm1)

! asize is max index for arrays y, yp, etc. nfill may not be an exact
! multiple of chunk, so the last block of emit, dryd etc could be smalle
! than chunk in size.
        asize = j1 - j0 + 1
!
! Zero temporary arrays
        em = 0.0
        dd = 0.0
        dw = 0.0
        m  = 0.0
        tc = 0.0
        ql = 0.0
        cs = 0.0
        ll = 0.0
        kj = 0.0   ! Set to minimum rate

! Assign appropriate parts of main arrays to temporary arrays
        DO j = 1, nc
          DO i = 1, asize
            k = j0 - 1 + i
            em(i,j) = emit(j,k)
            dd(i,j) = dryd(j,k)
            dw(i,j) = wetd(j,k)
          END DO
        END DO
        DO j = 1, ndj
          DO i = 1, asize
            k = j0 - 1 + i
            kj(i,j) = dj(j,k)
          END DO
        END DO
        DO i = 1, asize
          k = j0 - 1 + i
          m(i) = md(k)
          tc(i) = tl(k)
          ql(i) = q(k)
          ll(i) = liq(k)
          em(i,nc+1) = emit(nc+1,k)
          em(i,nc+2) = emit(nc+2,k)
        END DO
!
! If last part of arrays is < chunk, copy over parts to make the
! arrays m, tc and ql  full - avoids problems in CHEMCO and below with
! divide by zero errors.
        IF (asize < chunk) THEN
          m(asize+1:chunk) = md(1)
          tc(asize+1:chunk) = tl(1)
          ql(asize+1:chunk) = q(1)
        END IF

! Set cloud liquid water content, convert to vol water/vol air
        DO j = 1, asize
          ll(j) = ll(j) * m(j) * liq_const * one_over_cvf
        END DO

! Pick up the concentrations from the previous chemistry
! Fill in missing values if needed.
        DO i = 1, nc
          yp(1:asize,i) = xx(i,j0:j1) * m(1:asize)
          em(1:asize,i) = em(1:asize,i) * m(1:asize)
        END DO
        DO j = 1, asize
          em(j,nc+1) = em(j,nc+1) * m(1)
          em(j,nc+2) = em(j,nc+2) * m(1)
        END DO
        IF (asize < chunk) THEN
          DO i = 1, nc
            yp(asize+1:chunk,i) = yp(1,i)
            em(asize+1:chunk,i) = em(1,i)
          END DO
        END IF

! Set water vapour in molecules cm-3.
        DO j = 1, asize
          yp(j,i_h2o) = ql(j) * m(j)
        END DO

! Set concentrations to previous values for initialisation
        DO i = 1, nc
          DO j = 1, asize
            y(j,i) = yp(j,i)
          END DO
        END DO

! Find cells where the cloud liquid water content is greater than l_crit
        ngtlcrit = 0
        indx = 0
        DO j = 1, asize
          IF (ll(j) > l_crit) THEN
            ngtlcrit = ngtlcrit + 1
            indx(ngtlcrit) = j
          END IF
        END DO

! Set in-cloud concentrations if necessary, factor molec2mol converts fr
! molecules / cm3 to mol / l (air).

        DO i = 1, ngtlcrit
          j = indx(i)
          cth2o2(j) = y(j,i_h2o2) * molec2mol
          cthno3(j) = y(j,i_hno3) * molec2mol
          ctso2(j) = y(j,i_so2)   * molec2mol
          ctnh4(j) = y(j,i_ammsul)* molec2mol
          ctnh3(j) = y(j,i_nh3)   * molec2mol
          cs(j,10) = y(j,i_sa)    * molec2mol / ll(j) ! mol/l water
! Record the initial cloud concentrations.
          c1h2o2(j) = cth2o2(j)
          c1so2(j) = ctso2(j)
          c1so4(j) = cs(j,10)
          c1nh4(j) = ctnh4(j)
          c1nh3(j) = ctnh3(j)
        END DO

! Avoid ozone becoming too low.
        DO j = 1, asize
          k = j + j0 - 1
          IF (xx(i_o3,k) < small_conc) THEN
            xx(i_o3,k) = small_conc
            yp(j,i_o3) = xx(i_o3,k) * m(j)
            y(j,i_o3) = yp(j,i_o3)
          END IF
        END DO

! Set up stratospheric flag, true if in stratosphere.
! Ozone emission is only set if in the stratosphere.
        DO j = 1, asize
          in_strat(j) = (em(j,i_o3) /= 0.0)
        END DO

!   ***************
!   INITIALIZATION:
!   ***************

! Lagrangian time local starting time
        ctime = time

! Tmax local stop-time for integration
        tmax = ctime + astep

! Check photolysis rates are not < 0
        kj = MAX(kj, 1.0e-20)

! Calculate Rate Coefficients rc (cm, molecule, s units)
! Calculate Henry's law and equilibrium coefficients

! DEPENDS ON: chemco
        CALL CHEMCO(rc,tc,m,y(:,i_h2o),asize)
! DEPENDS ON: eqmcon
        CALL EQMCON(asize,tc,ke,kh)
! ------------------------------------

!   **************************************
!   Chemical time integration loop starts:
!   **************************************

! DTS is current chemistry step size.
! Set DTS depending on CO emission.
!O        IF(EM(8)>1.0E-19) THEN
!O          DTS=100.0
!O        ELSE
!O          DTS=300.0
!O        ENDIF

! Chemistry time step is constant for now
        dts = 300.0

! reus converts molecules/cm3 to atm.
         DO j = 1, asize
           reus(j) = molec2atm * tc(j)
         END DO

! Set all photolysis and reaction rate constants to zero
! if in the stratosphere.
         DO j = 1, asize
           IF (in_strat(j)) THEN
             kj(j,:) = 0.0
             rc(j,:) = 0.0
           END IF
         END DO

! Start calculation of new concentrations
        DO
          IF (ctime >= tmax) EXIT    ! Timestep loop start.
          ctime = ctime + dts
!
! Calculate new concentrations
!
! Set aqueous phase concentrations
          DO i = 1, ngtlcrit
            j = indx(i)
            cpso4(j) = cs(j,10)
            cph2o2(j) = cth2o2(j)
            cpso2(j) = ctso2(j)
            cpnh3(j) = ctnh3(j)
            cpnh4(j) = ctnh4(j)
            so4(j) = cpso4(j)
          END DO
          IF (ANY(cph2o2 < 0.0)) THEN
            WRITE(6,*) ' *** CPH2O2 -ve !, CTH2O2: ',cth2o2,            &
     &        ' CPH2O2: ',cph2o2
            WRITE(6,*) ' LAQ(3): ',laq(:,3)
            WRITE(6,*) ' HP: ',hp
            WRITE(6,*) ' CS: ',(cs(:,i),i=1,12)
          END IF
          IF (ANY(so4 < 0.0)) THEN
            WRITE(6,*) ' SO4 is negative, SO4: ',so4
          END IF

! Calculate the pH of the aqueous solution
! DEPENDS ON: phsol
          CALL PHSOL(asize,kh,ke,ll,tc,m,cthno3,ctso2,ctnh3,so4,hp,     &
     &      l_crit)

          f = 1.0
          DO i = 1, ngtlcrit
            j = indx(i)
            fso2(j) = (1.0 + ke(j,2) / hp(j))
            ph(j) = -LOG10(hp(j))
!      y(i_ph)=ph(j)*m(j) ! Multiply PH by M so that XX(74) comes out ri
          END DO

! Loop over backward-Euler iterations.
          DO k = 1, nit

!        -----------------------
!        AQUEOUS PHASE REACTIONS
!        -----------------------

            DO j = 1, asize
              IF (ll(j) > l_crit) THEN               ! Cloud water prese
!
! 1) Obtain dissolved species concentrations in (mol/l):
!  H2O2, O3, HNO3, SO2, NH3 from equilibrium equations.
!
                cs(j,1) = kh(j,1) * y(j,i_o3) * reus(j)           ! O3 (
                cs(j,2) = cthno3(j) / ll(j)                       ! HNO3
                cs(j,3) = cth2o2(j) / (ll(j)+pstar /                    &
     &            (kh(j,3) * rgc * tc(j) * 1.0e3))                ! H2O2
                cs(j,4) = ctso2(j) / (ll(j)+pstar /                     &
     &            (kh(j,4) * rgc * tc(j) * 1.0e3 * fso2(j)))      ! SO2
                cs(j,5) = ctnh3(j) / ll(j)                        ! NH3
                cs(j,6) = cs(j,2)                                 ! NO3
                cs(j,7) = cs(j,4) / (1.0 + hp(j) / ke(j,2))       ! HSO3
                cs(j,8) = cs(j,7) * ke(j,3) / hp(j)               ! SO3
                cs(j,9) = cs(j,5) / (1.0+ke(j,6)/(hp(j)*ke(j,4))) ! NH4
                cs(j,11) = kh(j,6) * xx_co2 * m(j) * reus(j)      ! CO2
                cs(j,12) = cs(j,11) / (1.0+hp(j) / ke(j,5))       ! HCO3

! 2) Calculate the production of SO4 and the loss of H2O2, and pass
!  these rates to the main model.

!    Fluxes:  (for budget printout)
!    FAQ(1) - SO2->SA (includes FAQ(2))
!    FAQ(2) - H2O2+SO2->SA
!    FAQ(3) - 2NH4+SO4->(NH4)2SO4
!    FAQ(4) - HSO3+O3->SO4
!    FAQ(5) - SO3+O3->SO4

!    Loss rates:
!    LAQ(6) - loss of in-cloud total NH3,
!    LAQ(7) - loss of in-cloud SO4,
!    LAQ(8) - loss of gaseous NH3,
!    LAQ(8) - gain of in-cloud total (NH4)2SO4,
!    LAQ(9) - loss of SA.
!    LAQ(10) - gain of Ammonium Sulphate.

                IF (cpnh3(j) > 1.0e-20 .AND. cpso4(j) > 1.0e-20) THEN
!    &            .AND. y(j,i_sa) > small_conc*m(j)) THEN
                  laq(j,6) = 2.0 * MIN(cs(j,9),cs(j,10)) * ll(j) /      &
     &              (dts * ctnh3(j))
                  laq(j,7) = MIN(cs(j,9),cs(j,10)) / (dts * cs(j,10))
                  laq(j,8) = laq(j,6) * ctnh3(j) / 2.0
                  laq(j,9) = laq(j,7) * cs(j,10) * ll(j) * mol2molec *  &
     &              cvf /  y(j,i_sa)
                  laq(j,10)= laq(j,9) * y(j,i_sa)
                ELSE
                  laq(j,6) = 0.0
                  laq(j,7) = 0.0
                  laq(j,8) = 0.0
                  laq(j,9) = 0.0
                  laq(j,10) = 0.0
                END IF
!
!        SO4(Aq)         CS(10)
! 0.1 is an empirical constant in mol/l
! (Martin & Damschen (1981),Bower et al., (1991))
                p(j) = cs(j,7)*cs(j,3)*rc(j,260)*hp(j) / (hp(j)+0.1)    &
     &            + cs(j,7)*cs(j,1)*rc(j,261)+cs(j,8)*cs(j,1)*rc(j,262)

!               write(6,'(a6,2x,es11.4)') 'p(so4)',p(j)
!               write(6,'(a6,2x,es11.4)') 'll    ',ll(j)
!               write(6,'(a6,2x,es11.4)') 'ctso2 ',ctso2(j)
!               write(6,'(a6,2x,es11.4)') 'rc 260',rc(j,260)
!               write(6,'(a6,2x,es11.4)') 'rc 261',rc(j,261)
!               write(6,'(a6,2x,es11.4)') 'rc 262',rc(j,262)
!
! scale production rate of SO4 so we don't oxidise more than 10% of the
! gaseous SO2 per timestep
                ee = p(j) * ll(j) * dts / ctso2(j)
                IF (ee > 0.1) THEN
                  f(j) = 0.1 / ee
                  p(j) = p(j) * f(j)
                ELSE
                  f(j) = 1.0
                END IF

                l(j) = laq(j,7)
                cs(j,10) = (cpso4(j) + dts*p(j)) / (1.0 + dts*l(j))

!         IF (ANY(cs(10,:) <= 0.0)) THEN
!           WRITE(6,*) 'SO4(aq) -ve'
!           WRITE(6,*) 'LAQ(7): ',LAQ(7,:)
!           WRITE(6,*) 'P: ',P
!           WRITE(6,*) 'HP ',HP
!           WRITE(6,*) 'CS(1) ',CS(1,:)
!           WRITE(6,*) 'CS(3) ',CS(3,:)
!           WRITE(6,*) 'CS(7) ',CS(7,:)
!           WRITE(6,*) 'CS(8) ',CS(8,:)
!           WRITE(6,*) 'RC(260) ',RC(260,:)
!           WRITE(6,*) 'RC(261) ',RC(261,:)
!           WRITE(6,*) 'RC(262) ',RC(262,:)
!           WRITE(6,*) 'Resetting to SA : 4.0e-18'
!         END IF
                IF (cs(j,10) < 0.0) THEN
                  yp(j,i_sa) = small_conc * m(j)
                  y(j,i_sa) = yp(j,i_sa)
                  cs(j,10) = y(j,i_sa) * molec2mol / ll(j)
                  xx(i_sa,j+j0-1) = small_conc
                END IF
!
!    Loss rates:
!    LAQ(2) - gain of SA,
!    LAQ(3) - loss of in-cloud total H2O2,
!    LAQ(4) - loss of gaseous H2O2,
!    LAQ(5) - loss of in-cloud total SO2.

! convert mol/l water to molecules/cm3 air
                laq(j,2) = p(j)*ll(j)*mol2molec*cvf
                laq(j,3) = ll(j)*f(j)*(cs(j,7)*cs(j,3)*rc(j,260)*       &
     &            hp(j) / (hp(j)+0.1)) / cth2o2(j)                 ! in
! convert mol/l water to molecules/cm3 air
                laq(j,4) = laq(j,3)*cth2o2(j)*mol2molec*cvf /           &
     &            y(j,i_h2o2)
                laq(j,5) = ll(j) * p(j) / ctso2(j)
                faq(j,1) = laq(j,2)
                faq(j,2) = laq(j,3) * cth2o2(j) * mol2molec * cvf
                faq(j,3) = laq(j,10)
                faq(j,4) = cs(j,1)*cs(j,7)*rc(j,261)*ll(j)*mol2molec*   &
     &            cvf*f(j)
                faq(j,5) = cs(j,1)*cs(j,8)*rc(j,262)*ll(j)*mol2molec*   &
     &            cvf*f(j)

!        Total cloud peroxide  CTH2O2
                p(j) = 0.0
                l(j) = laq(j,3)
!               cth2o2 = (cph2o2+dts*p)/(1.0+dts*l)
                cth2o2(j) = cph2o2(j) / (1.0 + dts * l(j))

!        Total cloud SO2  CTSO2
                p(j) = 0.0
                l(j) = laq(j,5)
!               ctso2 = (cpso2+dts*p)/(1.0+dts*l)
                ctso2(j) = cpso2(j) / (1.0 + dts * l(j))

!        Total cloud NH3  CTNH3
                p(j) = 0.0
                l(j) = laq(j,6)
!               ctnh3 = (cpnh3+dts*p)/(1.0+dts*l)
                ctnh3(j) = cpnh3(j) / (1.0 + dts * l(j))

!        Total cloud NH42SO4  CTNH4
                p(j) = laq(j,8)
                l(j) = 0.0
!           ctnh4 = (cpnh4+dts*p)/(1.0+dts*l)
                ctnh4(j) = cpnh4(j) + dts * p(j)

! No cloud water, set rates to zero
! Have to do each row separately (?) for vectorisation
              ELSE
                laq(j,1) = 0.0
                laq(j,2) = 0.0
                laq(j,3) = 0.0
                laq(j,4) = 0.0
                laq(j,5) = 0.0
                laq(j,6) = 0.0
                laq(j,7) = 0.0
                laq(j,8) = 0.0
                laq(j,9) = 0.0
                laq(j,10) = 0.0
                faq(j,1) = 0.0
                faq(j,2) = 0.0
                faq(j,3) = 0.0
                faq(j,4) = 0.0
                faq(j,5) = 0.0
              END IF
            END DO
!
!      Gaseous chemistry 70 species mechanism:
!      This section written automatically by MECHNEW from the file newme
!      with 178 equations.  6/II/97
!
!      Hand edits:              Double loss from four O3 reactions delet
!                               NO3 & N2O5 solved simultaneously.
!                               DDs put in, DW(39), DW(70) added.
!
!      If in stratosphere set concentrations of O(1D), O(3P) and OH
!      to zero, otherwise do normal chemistry calculation.
!
            DO j = 1, asize
              IF (in_strat(j)) THEN
                y(j,i_od) = 0.0
                y(j,i_op) = 0.0
                y(j,i_oh) = 0.0
              ELSE
!          OD           Y( 1)
        pt =                                                            &
     &+(kj(j,2)  *y(j,11))
        lt = 0.0                                                        &
     &+(rc(j,7)  )      +(rc(j,8)  *y(j,73))
        y(j, 1) = pt/lt
!
!          OP           y(:, 2)
        pt =                                                            &
     &+(kj(j,3)  *y(j,5 ))      +(kj(j,14) *y(j,6 ))                    &
     &+(rc(j,7)  *y(j,1 ))      +(kj(j,1)  *y(j,11))
        lt = 0.0                                                        &
     &+(rc(j,1)  )      +(rc(j,5)  *y(j,4 ))                            &
     &+(rc(j,16) *y(j,5 ))
        y(j, 2) = pt/lt
!
!          OH           y(j, 3)
        pt =                                                            &
     &+(kj(j,16) *y(j,52))      +(kj(j,16) *y(j,53))                    &
     &+(kj(j,16) *y(j,33))      +(kj(j,16) *y(j,35))                    &
     &+(kj(j,16) *y(j,22))      +(kj(j,16) *y(j,34))                    &
     &+(kj(j,4)  *y(j,14))      +(kj(j,5)  *y(j,13))                    &
     &+(rc(j,247)*y(j,53)*y(j,3 ))+(kj(j,4)  *y(j,14))                  &
     &+(rc(j,129)*y(j,11)*y(j,50)*0.36)                                 &
     &+(rc(j,245)*y(j,52)*y(j,3 ))                                      &
     &+(rc(j,123)*y(j,11)*y(j,28)*0.28)                                 &
     &+(rc(j,128)*y(j,11)*y(j,48)*0.27)                                 &
     &+(rc(j,102)*y(j,3 )*y(j,33))+(rc(j,104)*y(j,3 )*y(j,35))          &
     &+(rc(j,45) *y(j,3 )*y(j,22))+(rc(j,100)*y(j,3 )*y(j,34))          &
     &+(rc(j,17) *y(j,4 )*y(j,16))+(rc(j,34) *y(j,6 )*y(j,16))          &
     &+(rc(j,8)*y(j,1 )*y(j,73)*2.00)+(rc(j,14) *y(j,16)*y(j,11))
        lt =                                                            &
     &+(rc(j,248)*y(j,42))+(rc(j,249)*y(j,44))                          &
     &+(rc(j,251)*y(j,48))+(rc(j,253)*y(j,50))                          &
     &+(rc(j,236)*y(j,29))+(rc(j,237)*y(j,69))                          &
     &+(rc(j,245)*y(j,52))+(rc(j,247)*y(j,53))                          &
     &+(rc(j,230)*y(j,29))+(rc(j,232)*y(j,46))                          &
     &+(rc(j,234)*y(j,54))+(rc(j,235)*y(j,54))                          &
     &+(rc(j,214)*y(j,60))+(rc(j,215)*y(j,60))                          &
     &+(rc(j,224)*y(j,62))+(rc(j,225)*y(j,61))                          &
     &+(rc(j,125)*y(j,28))+(rc(j,204)*y(j,55))                          &
     &+(rc(j,206)*y(j,56))+(rc(j,210)*y(j,57))                          &
     &+(rc(j,100)*y(j,34))+(rc(j,102)*y(j,33))                          &
     &+(rc(j,104)*y(j,35))+(rc(j,109)*y(j,27))                          &
     &+(rc(j,86) *y(j,25))+(rc(j,92) *y(j,31))                          &
     &+(rc(j,94) *y(j,37))+(rc(j,98) *y(j,21))                          &
     &+(rc(j,70) *y(j,8 ))+(rc(j,71) *y(j,17))                          &
     &+(rc(j,75) *y(j,19))+(rc(j,81) *y(j,23))                          &
     &+(rc(j,44) *y(j,35))+(rc(j,59) *y(j,9 ))                          &
     &+(rc(j,63) *y(j,36))+(rc(j,66) *y(j,10))                          &
     &+(rc(j,44) *y(j,22))+(rc(j,45) *y(j,22))                          &
     &+(rc(j,44) *y(j,34))+(rc(j,44) *y(j,33))                          &
     &+(rc(j,31) *y(j,14))+(rc(j,33) *y(j,12))                          &
     &+(rc(j,35) *y(j,13))+(rc(j,39) *y(j,26))                          &
     &+(rc(j,13) *y(j,11))+(rc(j,21) *y(j,5 ))                          &
     &+(rc(j,24) *y(j,65))+(rc(j,30) *y(j,16))                          &
     &+(rc(j,136)*y(j,72))
        y(j, 3) = pt/lt
!
              END IF
            END DO ! End of stratospheric modification
!
!          NO           y(j, 4)
      DO j = 1, asize
        p(j) = em(j, 4)                                                 &
     &+(kj(j,3)  *y(j,5 ))      +(kj(j,13) *y(j,6 ))                    &
     &+(rc(j,217)*y(j,66)*y(j,5 ))+(rc(j,219)*y(j,64)*y(j,5 ))          &
     &+(rc(j,16) *y(j,5 )*y(j,2 ))+(rc(j,19) *y(j,5 )*y(j,6 ))
        l(j) = 0.0                                                      &
     &+(rc(j,254)*y(j,51))                                              &
     &+(rc(j,126)*y(j,43))+(rc(j,226)*y(j,63))                          &
     &+(rc(j,233)*y(j,47))+(rc(j,252)*y(j,49))                          &
     &+(rc(j,93) *y(j,32))+(rc(j,95) *y(j,38))                          &
     &+(rc(j,105)*y(j,40))+(rc(j,110)*y(j,41))                          &
     &+(rc(j,60) *y(j,15))+(rc(j,72) *y(j,18))                          &
     &+(rc(j,79) *y(j,20))+(rc(j,83) *y(j,24))                          &
     &+(rc(j,5)  *y(j,2 ))+(rc(j,11) *y(j,11))                          &
     &+(rc(j,15) *y(j,6 ))+(rc(j,17) *y(j,16))
        y(j, 4) = (yp(j, 4)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          NO2          y(j, 5)
      DO j = 1, asize
        p(j) = em(j, 5)                                                 &
     &+(kj(j,15) *y(j,7 ))      +(kj(j,17) *y(j,21))                    &
     &+(kj(j,10) *y(j,65))      +(kj(j,14) *y(j,6 ))                    &
     &+(rc(j,254)*y(j,51)*y(j,4 ))+(kj(j,5)  *y(j,13))                  &
     &+(rc(j,238)*y(j,6 )*y(j,69)*2.00)                                 &
     &+(rc(j,252)*y(j,49)*y(j,4 ))                                      &
     &+(rc(j,233)*y(j,47)*y(j,4 ))+(rc(j,237)*y(j,3 )*y(j,69))          &
     &+(rc(j,210)*y(j,57)*y(j,3 ))+(rc(j,226)*y(j,63)*y(j,4 ))          &
     &+(rc(j,204)*y(j,55)*y(j,3 ))+(rc(j,206)*y(j,56)*y(j,3 ))          &
     &+(rc(j,110)*y(j,41)*y(j,4 ))+(rc(j,126)*y(j,43)*y(j,4 ))          &
     &+(rc(j,95) *y(j,38)*y(j,4 ))+(rc(j,105)*y(j,4 )*y(j,40))          &
     &+(rc(j,83) *y(j,24)*y(j,4 ))+(rc(j,93) *y(j,32)*y(j,4 ))          &
     &+(rc(j,78) *y(j,21))      +(rc(j,79) *y(j,20)*y(j,4 ))            &
     &+(rc(j,60) *y(j,4 )*y(j,15))+(rc(j,72) *y(j,18)*y(j,4 ))          &
     &+(rc(j,29) *y(j,7 ))      +(rc(j,34) *y(j,6 )*y(j,16))            &
     &+(rc(j,24) *y(j,65)*y(j,3 ))                                      &
     &+(rc(j,27) *y(j,6 )*y(j,6 )*2.00)                                 &
     &+(rc(j,19) *y(j,5 )*y(j,6 ))+(rc(j,23) *y(j,65))                  &
     &+(rc(j,15) *y(j,4 )*y(j,6 )*2.00)                                 &
     &+(rc(j,17) *y(j,4 )*y(j,16))                                      &
     &+(rc(j,5)  *y(j,2 )*y(j,4 ))+(rc(j,11) *y(j,4 )*y(j,11))
        l(j) = 0.0 + dd(j,5)                                            &
     &+(rc(j,217)*y(j,66))+(rc(j,219)*y(j,64))                          &
     &+(rc(j,231)*y(j,45))+(kj(j,3)  )                                  &
     &+(rc(j,21) *y(j,3 ))+(rc(j,22) *y(j,16))                          &
     &+(rc(j,77) *y(j,20))+(rc(j,213)*y(j,59))                          &
     &+(rc(j,12) *y(j,11))+(rc(j,16) *y(j,2 ))                          &
     &+(rc(j,19) *y(j,6 ))+(rc(j,20) *y(j,6 ))
        y(j, 5) = (yp(j, 5)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          NO3          y(j, 6)
!      p(j) = em(j, 6)
!     &+(kj(j,15) *y(j,7 ))
!     &+(rc(j,35) *y(j,3 )*y(j,13))+(rc(j,98) *y(j,3 )*y(j,21))
!     &+(rc(j,12) *y(j,5 )*y(j,11))+(rc(j,29) *y(j,7 ))
!      l(j) = 0.0
!     &+(kj(j,13) )      +(kj(j,14) )
!     &+(rc(j,208)*y(j,19))+(rc(j,209)*y(j,48))
!     &+(rc(j,228)*y(j,60))+(rc(j,238)*y(j,69))
!     &+(rc(j,201)*y(j,17))+(rc(j,202)*y(j,23))
!     &+(rc(j,203)*y(j,27))+(rc(j,205)*y(j,28))
!     &+(rc(j,27) *y(j,6 ))+(rc(j,32) *y(j,16))
!     &+(rc(j,34) *y(j,16))+(rc(j,67) *y(j,10))
!     &+(rc(j,15) *y(j,4 ))+(rc(j,19) *y(j,5 ))
!     &+(rc(j,20) *y(j,5 ))+(rc(j,27) *y(j,6 ))
!      y(j, 6) = (yp(j, 6)+dts*p(j))/(1.0+dts*l(j))
!
!          N2O5         y(j, 7)
!      p(j) = em(j, 7)
!     &+(rc(j,20) *y(j,5 )*y(j,6 ))
!      l(j) = 0.0
!     &+(rc(j,29) )+(rc(j,43) )+(kj(j,15) )+(dw(j,7)  )
!      y(j, 7) = (yp(j, 7)+dts*p(j))/(1.0+dts*l(j))
!
!          NO3          y(j, 6) &          N2O5         y(j, 7)
      DO j = 1, asize
        p(j) =                                                          &
     &(rc(j,35) *y(j,3 )*y(j,13))                                       &
     &+(rc(j,98) *y(j,3 )*y(j,21))+(rc(j,12) *y(j,5 )*y(j,11))
        l(j) =                                                          &
     &(kj(j,13) )      +(kj(j,14) )                                     &
     &+(rc(j,208)*y(j,19))+(rc(j,209)*y(j,48))                          &
     &+(rc(j,228)*y(j,60))+(rc(j,238)*y(j,69))                          &
     &+(rc(j,201)*y(j,17))+(rc(j,202)*y(j,23))                          &
     &+(rc(j,203)*y(j,27))+(rc(j,205)*y(j,28))                          &
     &+(rc(j,27) *y(j,6 ))+(rc(j,32) *y(j,16))                          &
     &+(rc(j,34) *y(j,16))+(rc(j,67) *y(j,10))                          &
     &+(rc(j,15) *y(j,4 ))+(rc(j,19) *y(j,5 ))                          &
     &+(rc(j,20) *y(j,5 ))+(rc(j,27) *y(j,6 ))
        r1(j) = kj(j,15)+rc(j,29)
        r2(j) = (rc(j,20) *y(j,5 ))
        l1(j) = dd(j, 7)+dw(j, 7)                                       &
     &+(rc(j,29) )      +(rc(j,43) )      +(kj(j,15) )
        l2(j) = 1.0+l(j)*dts
        l3(j) = 1.0+l1(j)*dts
        y(j,6) = (l3(j)*(yp(j,6)+p(j)*dts)+r1(j)*dts*yp(j,7))/          &
     &    ((l3(j)*l2(j))-r1(j)*r2(j)*dts**2.0)
        y(j,7) = (yp(j,7)+r2(j)*dts*y(j,6))/l3(j)
      END DO
!
!          CO - total    y(j, 8)
      DO j = 1, asize
        p(j) = em(j, 8)                                                 &
     &+(kj(j,12) *y(j,44))                                              &
     &+(kj(j,8)  *y(j,19))      +(kj(j,11) *y(j,42))                    &
     &+(kj(j,6)  *y(j,10))      +(kj(j,7)  *y(j,10))                    &
     &+(rc(j,248)*y(j,42)*y(j,3 ))                                      &
     &+(rc(j,249)*y(j,44)*y(j,3 )*2.00)                                 &
     &+(rc(j,129)*y(j,11)*y(j,50)*0.76)                                 &
     &+(rc(j,223)*y(j,67)*y(j,10))                                      &
     &+(rc(j,124)*y(j,11)*y(j,28)*0.58)                                 &
     &+(rc(j,128)*y(j,11)*y(j,48)*0.78)                                 &
     &+(rc(j,112)*y(j,11)*y(j,27)*0.31)                                 &
     &+(rc(j,123)*y(j,11)*y(j,28)*0.40)                                 &
     &+(rc(j,66) *y(j,3 )*y(j,10))+(rc(j,67) *y(j,6 )*y(j,10))
        l(j) = 0.0 + dd(j,8)                                            &
     &+(rc(j,70) *y(j,3 ))
        y(j, 8) = (yp(j, 8)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!           CO(E) - emitted y(:,72)
      DO j = 1, asize
        p(j) = em(j,8)
        l(j) = 0.0 + dd(j,8)                                            &
     &+(rc(j,70) *y(j,3 ))
        y(j,72) = (yp(j,72)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH4          y(j, 9)
      DO j = 1, asize
        p(j) = em(j, 9)                                                 &
     &+(rc(j,123)*y(j,11)*y(j,28)*0.30)
        l(j) = 0.0                                                      &
     &+(rc(j,59) *y(j,3 ))
        y(j, 9) = (yp(j, 9)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          HCHO         y(j,10)
      DO j = 1, asize
        p(j) = em(j,10)                                                 &
     &+(kj(j,16) *y(j,52))      +(kj(j,16) *y(j,53))                    &
     &+(kj(j,12) *y(j,44))      +(kj(j,16) *y(j,22))                    &
     &+(rc(j,252)*y(j,49)*y(j,4 ))+(rc(j,254)*y(j,51)*y(j,4 ))          &
     &+(rc(j,245)*y(j,52)*y(j,3 ))+(rc(j,247)*y(j,53)*y(j,3 ))          &
!    &+(rc(j,242)*y(j,49)*y(j,15))+(rc(j,243)*y(j,51)*y(j,15))
     &+(rc(j,228)*y(j,6 )*y(j,60))+(rc(j,241)*y(j,47)*y(j,15))          &
     &+(rc(j,226)*y(j,63)*y(j,4 ))                                      &
     &+(rc(j,227)*y(j,63)*y(j,15)*2.00)                                 &
     &+(rc(j,204)*y(j,55)*y(j,3 ))+(rc(j,214)*y(j,3 )*y(j,60))          &
     &+(rc(j,128)*y(j,11)*y(j,48)*0.22)                                 &
     &+(rc(j,129)*y(j,11)*y(j,50)*0.24)                                 &
     &+(rc(j,126)*y(j,43)*y(j,4 ))                                      &
     &+(rc(j,127)*y(j,15)*y(j,43)*2.00)                                 &
     &+(rc(j,112)*y(j,11)*y(j,27)*0.47)                                 &
     &+(rc(j,123)*y(j,11)*y(j,28))                                      &
     &+(rc(j,111)*y(j,15)*y(j,41)*3.00)                                 &
     &+(rc(j,112)*y(j,11)*y(j,27))                                      &
     &+(rc(j,106)*y(j,15)*y(j,40))                                      &
     &+(rc(j,110)*y(j,41)*y(j,4 )*2.00)                                 &
     &+(rc(j,97) *y(j,32)*y(j,15))+(rc(j,98) *y(j,3 )*y(j,21))          &
     &+(rc(j,95) *y(j,38)*y(j,4 ))+(rc(j,96) *y(j,38)*y(j,15))          &
     &+(rc(j,80) *y(j,15)*y(j,20))+(rc(j,84) *y(j,24)*y(j,15))          &
     &+(rc(j,73) *y(j,18)*y(j,15))                                      &
     &+(rc(j,74) *y(j,15)*y(j,20)*2.00)                                 &
     &+(rc(j,62) *y(j,15)*y(j,15))+(rc(j,63) *y(j,36)*y(j,3 ))          &
     &+(rc(j,60) *y(j,4 )*y(j,15))                                      &
     &+(rc(j,61) *y(j,15)*y(j,15)*2.00)                                 &
!     &+(rc(j,40) *y(j,26)*y(j,15)) ! Don't think we want SO2+CH3O2
     &+(rc(j,45) *y(j,3 )*y(j,22))
        l(j) = 0.0 + dd(j,10)                                           &
     &+(kj(j,7)  )      +(dw(j,10) )                                    &
     &+(rc(j,66) *y(j,3 ))+(rc(j,67) *y(j,6 ))                          &
     &+(rc(j,223)*y(j,67))+(kj(j,6)  )
        y(j,10) = (yp(j,10)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          O3  total     y(j,11)
      DO j = 1, asize
        p(j) = em(j,11)                                                 &
     &+(rc(j,1)  *y(j,2 ))
        l(j) = 0.0 + dd(j,11)                                           &
     &+(kj(j,2)  )                                                      &
     &+(rc(j,129)*y(j,50))+(rc(j,216)*y(j,66))                          &
     &+(rc(j,218)*y(j,64))+(kj(j,1)  )                                  &
     &+(rc(j,124)*y(j,28))+(rc(j,128)*y(j,48))                          &
     &+(rc(j,112)*y(j,27))+(rc(j,123)*y(j,28))                          &
     &+(rc(j,11) *y(j,4 ))+(rc(j,12) *y(j,5 ))                          &
     &+(rc(j,13) *y(j,3 ))+(rc(j,14) *y(j,16))
        y(j,11) = (yp(j,11)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          O3(S) -strat y(j,71)
      DO j = 1, asize
        p(j) = em(j,11)
! Sinks due to losses of NO2, NO3 not accounted for.
        l(j) = 0.0 + dd(j,11)                                           &
     &+(rc(j,8)  *y(j,73)*y(j,1)/y(j,11))                               &
                                          ! O1D+H2O
     &+(rc(j,129)*y(j,50))+(rc(j,216)*y(j,66))                          &
     &+(rc(j,218)*y(j,64))                                              &
     &+(rc(j,124)*y(j,28))+(rc(j,128)*y(j,48))                          &
     &+(rc(j,112)*y(j,27))+(rc(j,123)*y(j,28))                          &
     &+(rc(j,13) *y(j,3 ))+(rc(j,14) *y(j,16))
        y(j,71) = (yp(j,71)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          H2           y(j,12)
      DO j = 1, asize
        p(j) = em(j,12)                                                 &
     &+(kj(j,7)  *y(j,10))                                              &
     &+(rc(j,112)*y(j,11)*y(j,27)*0.13)                                 &
     &+(rc(j,124)*y(j,11)*y(j,28)*0.24)
        l(j) = 0.0 + dd(j,12)                                           &
     &+(rc(j,33) *y(j,3 ))
        y(j,12) = (yp(j,12)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          HNO3         y(j,13)
      DO j = 1, asize
        p(j) = em(j,13)                                                 &
     &+(rc(j,228)*y(j,6 )*y(j,60))                                      &
     &+(rc(j,202)*y(j,6 )*y(j,23))+(rc(j,208)*y(j,6 )*y(j,19))          &
     &+(rc(j,67) *y(j,6 )*y(j,10))+(rc(j,201)*y(j,6 )*y(j,17))          &
     &+(rc(j,21) *y(j,5 )*y(j,3 ))+(rc(j,32) *y(j,6 )*y(j,16))
        l(j) = 0.0 + dd(j,13)                                           &
     &+(rc(j,35) *y(j,3 ))+(rc(j,43) ) +(kj(j,5)  ) +(dw(j,13) )
        y(j,13) = (yp(j,13)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          H2O2         y(j,14)
      DO j = 1, asize
        p(j) = em(j,14)                                                 &
     &+(rc(j,36) *y(j,16)*y(j,16))
        l(j) = 0.0 + dd(j,14) + laq(j,4)                                &
     &+(rc(j,31) *y(j,3))+(kj(j,4))      +(dw(j,14) )
        y(j,14) = (yp(j,14)+dts*p(j))/(1.0+dts*l(j))
!       y(j,14) = MAX(y(j,14),1.0) ! stop -ve vals; don't need this.
      END DO
!
!          CH3O2        y(j,15)
      DO j = 1, asize
        p(j) = em(j,15)                                                 &
     &+(kj(j,8)  *y(j,19))      +(kj(j,9)  *y(j,37))                    &
     &+(rc(j,220)*y(j,64))      +(rc(j,222)*y(j,67))                    &
     &+(rc(j,91) *y(j,20)*y(j,20)*2.00)                                 &
     &+(rc(j,123)*y(j,11)*y(j,28)*0.58)                                 &
     &+(rc(j,79) *y(j,20)*y(j,4 ))+(rc(j,80) *y(j,15)*y(j,20))          &
     &+(rc(j,44) *y(j,3 )*y(j,22))+(rc(j,59) *y(j,3 )*y(j,9 ))
        l(j) = 0.0                                                      &
     &+(rc(j,227)*y(j,63))+(rc(j,241)*y(j,47))                          &
!    &+(rc(j,242)*y(j,49))+(rc(j,243)*y(j,51))
     &+(rc(j,97) *y(j,32))+(rc(j,106)*y(j,40))                          &
     &+(rc(j,111)*y(j,41))+(rc(j,127)*y(j,43))                          &
     &+(rc(j,74) *y(j,20))+(rc(j,80) *y(j,20))                          &
     &+(rc(j,84) *y(j,24))+(rc(j,96) *y(j,38))                          &
     &+(rc(j,62) *y(j,15))+(rc(j,62) *y(j,15))                          &
     &+(rc(j,65) *y(j,16))+(rc(j,73) *y(j,18))                          &
     &+(rc(j,60) *y(j,4 ))                                              &
     &+(rc(j,61) *y(j,15))+(rc(j,61) *y(j,15))
        y(j,15) = (yp(j,15)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          HO2          y(j,16)
      DO j = 1, asize
        p(j) = em(j,16)                                                 &
     &+(kj(j,16) *y(j,52))      +(kj(j,16) *y(j,53))                    &
     &+(kj(j,16) *y(j,33))      +(kj(j,16) *y(j,35))                    &
     &+(kj(j,16) *y(j,22))      +(kj(j,16) *y(j,34))                    &
     &+(kj(j,10) *y(j,65))      +(kj(j,11) *y(j,42))                    &
     &+(kj(j,6)  *y(j,10))      +(kj(j,8)  *y(j,19))                    &
     &+(rc(j,254)*y(j,51)*y(j,4 ))+(kj(j,6)  *y(j,10))                  &
     &+(rc(j,249)*y(j,44)*y(j,3 ))+(rc(j,252)*y(j,49)*y(j,4 ))          &
!    &+(rc(j,242)*y(j,49)*y(j,15)*2.00)
!    &+(rc(j,243)*y(j,51)*y(j,15)*2.00)
     &+(rc(j,234)*y(j,3 )*y(j,54))                                      &
     &+(rc(j,241)*y(j,47)*y(j,15)*2.00)                                 &
     &+(rc(j,230)*y(j,3 )*y(j,29))+(rc(j,233)*y(j,47)*y(j,4 ))          &
     &+(rc(j,224)*y(j,3 )*y(j,62))+(rc(j,227)*y(j,63)*y(j,15))          &
     &+(rc(j,215)*y(j,3 )*y(j,60))+(rc(j,223)*y(j,67)*y(j,10))          &
     &+(rc(j,205)*y(j,6 )*y(j,28))+(rc(j,209)*y(j,6 )*y(j,48))          &
     &+(rc(j,129)*y(j,11)*y(j,50)*0.36)                                 &
     &+(rc(j,203)*y(j,6 )*y(j,27))                                      &
     &+(rc(j,127)*y(j,15)*y(j,43)*2.00)                                 &
     &+(rc(j,128)*y(j,11)*y(j,48)*0.27)                                 &
     &+(rc(j,124)*y(j,11)*y(j,28)*0.18)                                 &
     &+(rc(j,126)*y(j,43)*y(j,4 ))                                      &
     &+(rc(j,112)*y(j,11)*y(j,27)*0.20)                                 &
     &+(rc(j,123)*y(j,11)*y(j,28)*0.30)                                 &
     &+(rc(j,110)*y(j,41)*y(j,4 ))                                      &
     &+(rc(j,111)*y(j,15)*y(j,41)*2.00)                                 &
     &+(rc(j,97) *y(j,32)*y(j,15))+(rc(j,106)*y(j,15)*y(j,40))          &
     &+(rc(j,93) *y(j,32)*y(j,4 ))+(rc(j,96) *y(j,38)*y(j,15))          &
     &+(rc(j,84) *y(j,24)*y(j,15)*2.00)                                 &
     &+(rc(j,90) *y(j,18)*y(j,18)*2.00)                                 &
     &+(rc(j,80) *y(j,15)*y(j,20))+(rc(j,83) *y(j,24)*y(j,4 ))          &
     &+(rc(j,72) *y(j,18)*y(j,4 ))                                      &
     &+(rc(j,73) *y(j,18)*y(j,15)*2.00)                                 &
     &+(rc(j,67) *y(j,6 )*y(j,10))+(rc(j,70) *y(j,3 )*y(j,8 ))          &
     &+(rc(j,63) *y(j,36)*y(j,3 ))+(rc(j,66) *y(j,3 )*y(j,10))          &
     &+(rc(j,60) *y(j,4 )*y(j,15))                                      &
     &+(rc(j,61) *y(j,15)*y(j,15)*2.00)                                 &
     &+(rc(j,39) *y(j,3 )*y(j,26))                                      &
     &+(rc(j,31) *y(j,3 )*y(j,14))+(rc(j,33) *y(j,3 )*y(j,12))          &
     &+(rc(j,13) *y(j,3 )*y(j,11))+(rc(j,23) *y(j,65))
        l(j) = 0.0                                                      &
     &+(rc(j,246)*y(j,51))                                              &
     &+(rc(j,212)*y(j,59))+(rc(j,221)*y(j,67))                          &
     &+(rc(j,240)*y(j,45))+(rc(j,244)*y(j,49))                          &
     &+(rc(j,65) *y(j,15))+(rc(j,99) *y(j,18))                          &
     &+(rc(j,101)*y(j,32))+(rc(j,103)*y(j,24))                          &
     &+(rc(j,32) *y(j,6 ))+(rc(j,34) *y(j,6 ))                          &
     &+(rc(j,36) *y(j,16))+(rc(j,36) *y(j,16))                          &
     &+(rc(j,14) *y(j,11))+(rc(j,17) *y(j,4 ))                          &
     &+(rc(j,22) *y(j,5 ))+(rc(j,30) *y(j,3 ))
        y(j,16) = (yp(j,16)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          C2H6         y(j,17)
      DO j = 1, asize
        p(j) = em(j,17)
        l(j) = 0.0                                                      &
     &+(rc(j,71) *y(j,3 ))+(rc(j,201)*y(j,6 ))
        y(j,17) = (yp(j,17)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          C2H5O2       y(j,18)
      DO j = 1, asize
        p(j) = em(j,18)                                                 &
     &+(rc(j,201)*y(j,6 )*y(j,17))+(kj(j,9)  *y(j,25))                  &
     &+(rc(j,44) *y(j,3 )*y(j,34))+(rc(j,71) *y(j,3 )*y(j,17))
        l(j) = 0.0                                                      &
     &+(rc(j,99) *y(j,16))                                              &
     &+(rc(j,72) *y(j,4 ))+(rc(j,73) *y(j,15))                          &
     &+(rc(j,90) *y(j,18))+(rc(j,90) *y(j,18))
        y(j,18) = (yp(j,18)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH3CHO       y(j,19)
      DO j = 1, asize
        p(j) = em(j,19)                                                 &
     &+(kj(j,16) *y(j,34))                                              &
     &+(rc(j,127)*y(j,15)*y(j,43))+(rc(j,206)*y(j,56)*y(j,3 ))          &
     &+(rc(j,124)*y(j,11)*y(j,28))+(rc(j,126)*y(j,43)*y(j,4 ))          &
     &+(rc(j,105)*y(j,4 )*y(j,40))+(rc(j,106)*y(j,15)*y(j,40))          &
     &+(rc(j,90) *y(j,18)*y(j,18)*2.00)                                 &
     &+(rc(j,100)*y(j,3 )*y(j,34))                                      &
     &+(rc(j,72) *y(j,18)*y(j,4 ))+(rc(j,73) *y(j,18)*y(j,15))
        l(j) = 0.0 + dd(j,19)                                           &
     &+(rc(j,75) *y(j,3 ))+(rc(j,208)*y(j,6 ))+(kj(j,8)  )
        y(j,19) = (yp(j,19)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH3COO2      y(j,20)
      DO j = 1, asize
        p(j) = em(j,20)                                                 &
     &+(kj(j,11) *y(j,42))      +(kj(j,17) *y(j,21))                    &
     &+(kj(j,9)  *y(j,25))      +(kj(j,9)  *y(j,37))                    &
     &+(rc(j,208)*y(j,6 )*y(j,19))+(rc(j,248)*y(j,42)*y(j,3 ))          &
     &+(rc(j,105)*y(j,4 )*y(j,40))+(rc(j,106)*y(j,15)*y(j,40))          &
     &+(rc(j,95) *y(j,38)*y(j,4 ))+(rc(j,96) *y(j,38)*y(j,15))          &
     &+(rc(j,75) *y(j,3 )*y(j,19))+(rc(j,78) *y(j,21))
        l(j) = 0.0                                                      &
     &+(rc(j,91) *y(j,20))+(rc(j,91) *y(j,20))                          &
     &+(rc(j,74) *y(j,15))+(rc(j,77) *y(j,5 ))                          &
     &+(rc(j,79) *y(j,4 ))+(rc(j,80) *y(j,15))
        y(j,20) = (yp(j,20)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          PAN          y(j,21)
      DO j = 1, asize
        p(j) = em(j,21)                                                 &
     &+(rc(j,77) *y(j,20)*y(j,5 ))
        l(j) = 0.0 + dd(j,21)                                           &
     &+(rc(j,78) )      +(rc(j,98) *y(j,3 ))+(kj(j,17) )
        y(j,21) = (yp(j,21)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH3OOH       y(j,22)
      DO j = 1, asize
        p(j) = em(j,22)                                                 &
     &+(rc(j,65) *y(j,15)*y(j,16))
        l(j) = 0.0 + dd(j,22)                                           &
     &+(rc(j,44) *y(j,3 ))+(rc(j,45) *y(j,3 ))                          &
     &+(kj(j,16) )      +(dw(j,22) )
        y(j,22) = (yp(j,22)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          NC4H10       y(j,23)
      DO j = 1, asize
        p(j) = em(j,23)
        l(j) = 0.0                                                      &
     &+(rc(j,81) *y(j,3 ))+(rc(j,202)*y(j,6 ))
        y(j,23) = (yp(j,23)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          SC4H9O2      y(j,24)
      DO j = 1, asize
        p(j) = em(j,24)                                                 &
     &+(rc(j,202)*y(j,6 )*y(j,23))                                      &
     &+(rc(j,44) *y(j,3 )*y(j,35))                                      &
     &+(rc(j,81) *y(j,3 )*y(j,23))
        l(j) = 0.0                                                      &
     &+(rc(j,83) *y(j,4 ))+(rc(j,84) *y(j,15))                          &
     &+(rc(j,103)*y(j,16))
        y(j,24) = (yp(j,24)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH3COE       y(j,25)
      DO j = 1, asize
        p(j) = em(j,25)                                                 &
     &+(rc(j,104)*y(j,3 )*y(j,35))+(kj(j,16) *y(j,35))                  &
     &+(rc(j,83) *y(j,24)*y(j,4 ))+(rc(j,84) *y(j,24)*y(j,15))
        l(j) = 0.0                                                      &
     &+(rc(j,86) *y(j,3 ))+(kj(j,9)  )
        y(j,25) = (yp(j,25)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          SO2          y(j,26)
      DO j = 1, asize
        p(j) = em(j,26)                                                 &
     &+(rc(j,220)*y(j,64))
        l(j) = 0.0 + dd(j,26)                                           &
     &+(rc(j,39) *y(j,3 ))                                              &
     &+(dw(j,26) )
        y(j,26) = (yp(j,26)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          C2H4         y(j,27)
      DO j = 1, asize
        p(j) = em(j,27)
        l(j) = 0.0                                                      &
     &+(rc(j,109)*y(j,3 ))+(rc(j,112)*y(j,11))                          &
     &+(rc(j,203)*y(j,6 ))
        y(j,27) = (yp(j,27)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          C3H6         y(j,28)
      DO j = 1, asize
        p(j) = em(j,28)
        l(j) = 0.0                                                      &
     &+(rc(j,205)*y(j,6 ))                                              &
     &+(rc(j,123)*y(j,11))+(rc(j,124)*y(j,11))                          &
     &+(rc(j,125)*y(j,3 ))
        y(j,28) = (yp(j,28)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!
!          OXYL         y(j,29)
      DO j = 1, asize
        p(j) = em(j,29)
        l(j) = 0.0                                                      &
     &+(rc(j,230)*y(j,3 ))+(rc(j,236)*y(j,3 ))
        y(j,29) = (yp(j,29)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          SA           y(j,30)
      DO j = 1, asize
        p(j) = em(j,30)                                                 &
     &+(rc(j,222)*y(j,67))+laq(j,2)                                     &
     &+(rc(j,39) *y(j,3 )*y(j,26))
        l(j) = 0.0 + dd(j,30)                                           &
     &+(dw(j,30) ) + laq(j,9)
        y(j,30) = (yp(j,30)+dts*p(j))/(1.0+dts*l(j))
        y(j,30) = MAX(y(j,30),1.0)  ! stop -ve vals, div by zero
      END DO
!
!          C3H8         y(j,31)
      DO j = 1, asize
        p(j) = em(j,31)
        l(j) = 0.0                                                      &
     &+(rc(j,92) *y(j,3 ))
        y(j,31) = (yp(j,31)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          C3H7O2       y(j,32)
      DO j = 1, asize
        p(j) = em(j,32)                                                 &
     &+(rc(j,44) *y(j,3 )*y(j,33))+(rc(j,92) *y(j,3 )*y(j,31))
        l(j) = 0.0                                                      &
     &+(rc(j,93) *y(j,4 ))+(rc(j,97) *y(j,15))                          &
     &+(rc(j,101)*y(j,16))
        y(j,32) = (yp(j,32)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          C3H7OOH      y(j,33)
      DO j = 1, asize
        p(j) = em(j,33)                                                 &
     &+(rc(j,101)*y(j,16)*y(j,32))
        l(j) = 0.0 + dd(j,33)                                           &
     &+(rc(j,44) *y(j,3 ))+(rc(j,102)*y(j,3 ))                          &
     &+(kj(j,16) )      +(dw(j,33) )
        y(j,33) = (yp(j,33)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          C2H5OOH      y(j,34)
      DO j = 1, asize
        p(j) = em(j,34)                                                 &
     &+(rc(j,99) *y(j,16)*y(j,18))
        l(j) = 0.0 + dd(j,34)                                           &
     &+(rc(j,44) *y(j,3 ))+(rc(j,100)*y(j,3 ))                          &
     &+(kj(j,16) )      +(dw(j,34) )
        y(j,34) = (yp(j,34)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          C4H9OOH      y(j,35)
      DO j = 1, asize
        p(j) = em(j,35)                                                 &
     &+(rc(j,103)*y(j,16)*y(j,24))
        l(j) = 0.0 + dd(j,35)                                           &
     &+(rc(j,44) *y(j,3 ))+(rc(j,104)*y(j,3 ))                          &
     &+(kj(j,16) )      +(dw(j,35) )
        y(j,35) = (yp(j,35)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH3OH        y(j,36)
      DO j = 1, asize
        p(j) = em(j,36)                                                 &
     &+(rc(j,62) *y(j,15)*y(j,15))                                      &
     &+(rc(j,123)*y(j,11)*y(j,28)*0.12)
        l(j) = 0.0 + dd(j,36)                                           &
     &+(rc(j,63) *y(j,3 ))
        y(j,36) = (yp(j,36)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          ACETONE      y(j,37)
      DO j = 1, asize
        p(j) = em(j,37)                                                 &
     &+(rc(j,102)*y(j,3 )*y(j,33))+(kj(j,16) *y(j,33))                  &
     &+(rc(j,93) *y(j,32)*y(j,4 ))+(rc(j,97) *y(j,32)*y(j,15))
        l(j) = 0.0 + dd(j,37)                                           &
     &+(rc(j,94) *y(j,3 ))+(kj(j,9)  )
        y(j,37) = (yp(j,37)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          ACETO2       y(j,38)
      DO j = 1, asize
        p(j) = em(j,38)                                                 &
     &+(rc(j,94) *y(j,3 )*y(j,37))
        l(j) = 0.0                                                      &
     &+(rc(j,95) *y(j,4 ))+(rc(j,96) *y(j,15))
        y(j,38) = (yp(j,38)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          NH42SO4      y(j,39)
      DO j = 1, asize
        p(j) = em(j,39)
        l(j) = 0.0 + dd(j,39) + dw(j,39)
        y(j,39) = (yp(j,39)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH3COX       y(j,40)
      DO j = 1, asize
        p(j) = em(j,40)                                                 &
     &+(rc(j,86) *y(j,3 )*y(j,25))
        l(j) = 0.0                                                      &
     &+(rc(j,105)*y(j,4 ))+(rc(j,106)*y(j,15))
        y(j,40) = (yp(j,40)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH2O2C       y(j,41)
      DO j = 1, asize
        p(j) = em(j,41)                                                 &
     &+(rc(j,109)*y(j,3 )*y(j,27))
        l(j) = 0.0                                                      &
     &+(rc(j,110)*y(j,4 ))+(rc(j,111)*y(j,15))
        y(j,41) = (yp(j,41)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          MGLYOX       y(j,42)
      DO j = 1, asize
        p(j) = em(j,42)                                                 &
     &+(kj(j,16) *y(j,53))                                              &
     &+(rc(j,247)*y(j,53)*y(j,3 ))+(rc(j,254)*y(j,51)*y(j,4 ))          &
     &+(rc(j,241)*y(j,47)*y(j,15))                                      &
!    &+(rc(j,243)*y(j,51)*y(j,15))
     &+(rc(j,233)*y(j,47)*y(j,4 ))+(rc(j,240)*y(j,45)*y(j,16))          &
     &+(rc(j,129)*y(j,11)*y(j,50))+(rc(j,230)*y(j,3 )*y(j,29))
        l(j) = 0.0                                                      &
     &+(rc(j,248)*y(j,3 ))+(kj(j,11) )
        y(j,42) = (yp(j,42)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH3CHX       y(j,43)
      DO j = 1, asize
        p(j) = em(j,43)                                                 &
     &+(rc(j,125)*y(j,3 )*y(j,28))
        l(j) = 0.0                                                      &
     &+(rc(j,126)*y(j,4 ))+(rc(j,127)*y(j,15))
        y(j,43) = (yp(j,43)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          GLYOX        y(j,44)
      DO j = 1, asize
        p(j) = em(j,44)                                                 &
     &+(rc(j,238)*y(j,6 )*y(j,69))+(rc(j,241)*y(j,47)*y(j,15))          &
     &+(rc(j,234)*y(j,3 )*y(j,54))+(rc(j,237)*y(j,3 )*y(j,69))          &
     &+(rc(j,212)*y(j,16)*y(j,59))+(rc(j,233)*y(j,47)*y(j,4 ))
        l(j) = 0.0                                                      &
     &+(rc(j,249)*y(j,3 ))+(kj(j,12) )
        y(j,44) = (yp(j,44)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          OXYL1        y(j,45)
      DO j = 1, asize
        p(j) = em(j,45)                                                 &
     &+(rc(j,236)*y(j,3 )*y(j,29))
        l(j) = 0.0                                                      &
     &+(rc(j,231)*y(j,5 ))+(rc(j,240)*y(j,16))
        y(j,45) = (yp(j,45)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          MEMALD       y(j,46)
      DO j = 1, asize
        p(j) = em(j,46)                                                 &
     &+(rc(j,238)*y(j,6 )*y(j,69))+(rc(j,240)*y(j,45)*y(j,16))          &
     &+(rc(j,234)*y(j,3 )*y(j,54))+(rc(j,237)*y(j,3 )*y(j,69))          &
     &+(rc(j,212)*y(j,16)*y(j,59))+(rc(j,230)*y(j,3 )*y(j,29))
        l(j) = 0.0                                                      &
     &+(rc(j,232)*y(j,3 ))
        y(j,46) = (yp(j,46)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          MEMALD1      y(j,47)
      DO j = 1, asize
        p(j) = em(j,47)                                                 &
     &+(rc(j,232)*y(j,3 )*y(j,46))
        l(j) = 0.0                                                      &
     &+(rc(j,233)*y(j,4 ))+(rc(j,241)*y(j,15))
        y(j,47) = (yp(j,47)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          C5H8         y(j,48)
      DO j = 1, asize
        p(j) = em(j,48)
        l(j) = 0.0                                                      &
     &+(rc(j,128)*y(j,11))+(rc(j,209)*y(j,6 ))                          &
     &+(rc(j,251)*y(j,3 ))
        y(j,48) = (yp(j,48)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          RO2IP1       y(j,49)
      DO j = 1, asize
        p(j) = em(j,49)                                                 &
     &+(rc(j,251)*y(j,3 )*y(j,48))
        l(j) = 0.0                                                      &
!    &+(rc(j,242)*y(j,15))
     &+(rc(j,244)*y(j,16))                                              &
     &+(rc(j,252)*y(j,4 ))
        y(j,49) = (yp(j,49)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          MVK          y(j,50)
      DO j = 1, asize
        p(j) = em(j,50)                                                 &
     &+(rc(j,252)*y(j,49)*y(j,4 ))+(kj(j,16) *y(j,52))                  &
     &+(rc(j,242)*y(j,49)*y(j,15))                                      &
!    &+(rc(j,245)*y(j,52)*y(j,3 ))
     &+(rc(j,128)*y(j,11)*y(j,48))+(rc(j,210)*y(j,57)*y(j,3 ))
        l(j) = 0.0                                                      &
     &+(rc(j,129)*y(j,11))+(rc(j,253)*y(j,3 ))
        y(j,50) = (yp(j,50)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          RO2IP2       y(j,51)
      DO j = 1, asize
        p(j) = em(j,51)                                                 &
     &+(rc(j,253)*y(j,3 )*y(j,50))
        l(j) = 0.0                                                      &
!    &+(rc(j,243)*y(j,15))
     &+(rc(j,246)*y(j,16))                                              &
     &+(rc(j,254)*y(j,4 ))
        y(j,51) = (yp(j,51)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          ISOPOOH      y(j,52)
      DO j = 1, asize
        p(j) = em(j,52)                                                 &
     &+(rc(j,244)*y(j,49)*y(j,16))
        l(j) = 0.0 + dd(j,52)                                           &
     &+(rc(j,245)*y(j,3 ))+(kj(j,16) )      +(dw(j,52) )
        y(j,52) = (yp(j,52)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          MVKOOH       y(j,53)
      DO j = 1, asize
        p(j) = em(j,53)                                                 &
     &+(rc(j,246)*y(j,51)*y(j,16))
        l(j) = 0.0 + dd(j,53)                                           &
     &+(rc(j,247)*y(j,3 ))+(kj(j,16) )      +(dw(j,53) )
        y(j,53) = (yp(j,53)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          TOLUEN       y(j,54)
      DO j = 1, asize
        p(j) = em(j,54)
        l(j) = 0.0                                                      &
     &+(rc(j,234)*y(j,3 ))+(rc(j,235)*y(j,3 ))
        y(j,54) = (yp(j,54)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          RNC2H4       y(j,55)
      DO j = 1, asize
        p(j) = em(j,55)                                                 &
     &+(rc(j,203)*y(j,6 )*y(j,27))
        l(j) = 0.0                                                      &
     &+(rc(j,204)*y(j,3 ))
        y(j,55) = (yp(j,55)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          RNC3H6       y(j,56)
      DO j = 1, asize
        p(j) = em(j,56)                                                 &
     &+(rc(j,205)*y(j,6 )*y(j,28))
        l(j) = 0.0                                                      &
     &+(rc(j,206)*y(j,3 ))
        y(j,56) = (yp(j,56)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          RNC5H8       y(j,57)
      DO j = 1, asize
        p(j) = em(j,57)                                                 &
     &+(rc(j,209)*y(j,6 )*y(j,48))
        l(j) = 0.0                                                      &
     &+(rc(j,210)*y(j,3 ))
        y(j,57) = (yp(j,57)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          NAER         y(j,58)
      DO j = 1, asize
        p(j) = em(j,58)                                                 &
     &+(rc(j,43) *y(j,13))      +(rc(j,43) *y(j,69))                    &
     &+(rc(j,43) *y(j,7 ))      +(rc(j,43) *y(j,7 ))
        l(j) = 0.0 + dd(j,58)                                           &
     &+(dw(j,58) )
        y(j,58) = (yp(j,58)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          TOLP1        y(j,59)
      DO j = 1, asize
        p(j) = em(j,59)                                                 &
     &+(rc(j,235)*y(j,3 )*y(j,54))
        l(j) = 0.0                                                      &
     &+(rc(j,212)*y(j,16))+(rc(j,213)*y(j,5 ))
        y(j,59) = (yp(j,59)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          DMS          y(j,60)
      DO j = 1, asize
        p(j) = em(j,60)
        l(j) = 0.0 + dd(j,60)                                           &
     &+(rc(j,214)*y(j,3 ))+(rc(j,215)*y(j,3 ))                          &
     &+(rc(j,228)*y(j,6 ))
        y(j,60) = (yp(j,60)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          DMSO2        y(j,61)
      DO j = 1, asize
        p(j) = em(j,61)                                                 &
     &+(rc(j,224)*y(j,3 )*y(j,62))
        l(j) = 0.0                                                      &
     &+(rc(j,225)*y(j,3 ))
        y(j,61) = (yp(j,61)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          DMSO         y(j,62)
      DO j = 1, asize
        p(j) = em(j,62)                                                 &
     &+(rc(j,215)*y(j,3 )*y(j,60))
        l(j) = 0.0                                                      &
     &+(rc(j,224)*y(j,3 ))
        y(j,62) = (yp(j,62)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          DMSP         y(j,63)
      DO j = 1, asize
        p(j) = em(j,63)                                                 &
     &+(rc(j,225)*y(j,3 )*y(j,61))
        l(j) = 0.0                                                      &
     &+(rc(j,226)*y(j,4 ))+(rc(j,227)*y(j,15))
        y(j,63) = (yp(j,63)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH3SO2       y(j,64)
      DO j = 1, asize
        p(j) = em(j,64)                                                 &
     &+(rc(j,226)*y(j,63)*y(j,4 ))+(rc(j,227)*y(j,63)*y(j,15))          &
     &+(rc(j,216)*y(j,66)*y(j,11))+(rc(j,217)*y(j,66)*y(j,5 ))
        l(j) = 0.0                                                      &
     &+(rc(j,218)*y(j,11))+(rc(j,219)*y(j,5 ))+(rc(j,220))
        y(j,64) = (yp(j,64)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          HO2NO2       y(j,65)
      DO j = 1, asize
        p(j) = em(j,65)                                                 &
     &+(rc(j,22) *y(j,5 )*y(j,16))
        l(j) = 0.0                                                      &
     &+(rc(j,23) )      +(rc(j,24) *y(j,3 ))+(kj(j,10) )
        y(j,65) = (yp(j,65)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH3SO        y(j,66)
      DO j = 1, asize
        p(j) = em(j,66)                                                 &
     &+(rc(j,214)*y(j,3 )*y(j,60))+(rc(j,228)*y(j,6 )*y(j,60))
        l(j) = 0.0                                                      &
     &+(rc(j,216)*y(j,11))+(rc(j,217)*y(j,5 ))
        y(j,66) = (yp(j,66)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          CH3SO3       y(j,67)
      DO j = 1, asize
        p(j) = em(j,67)                                                 &
     &+(rc(j,218)*y(j,64)*y(j,11))+(rc(j,219)*y(j,64)*y(j,5 ))
        l(j) = 0.0                                                      &
     &+(rc(j,221)*y(j,16))+(rc(j,222))      +(rc(j,223)*y(j,10))
        y(j,67) = (yp(j,67)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          MSA          y(j,68)
      DO j = 1, asize
        p(j) = em(j,68)                                                 &
     &+(rc(j,221)*y(j,67)*y(j,16))+(rc(j,223)*y(j,67)*y(j,10))
        l(j) = 0.0 + dd(j,68)                                           &
     &+(dw(j,68) )
        y(j,68) = (yp(j,68)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          ORGNIT       y(j,69)
      DO j = 1, asize
        p(j) = em(j,69)                                                 &
     &+(rc(j,213)*y(j,5 )*y(j,59))+(rc(j,231)*y(j,45)*y(j,5 ))
        l(j) = 0.0 + dd(j,69)                                           &
     &+(rc(j,43) )      +(rc(j,237)*y(j,3 ))                            &
     &+(rc(j,238)*y(j,6 ))+(dw(j,69) )
        y(j,69) = (yp(j,69)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          NH3          y(j,70)
      DO j = 1, asize
        p(j) = em(j,70)
        l(j) = 0.0 + dd(j,70) + dw(j,70)
        y(j,70) = (yp(j,70)+dts*p(j))/(1.0+dts*l(j))
      END DO
!
!          BE7          y(j,74)
      DO j = 1, asize
        p(j) = em(j,I_BE7)
        l(j) = 0. + dd(j,I_BE7) + dw(j,I_BE7)                           &
     & +rc(j,290)
        y(j,I_BE7) = (yp(j,I_BE7)+dts*p(j))/(1.0+dts*l(j))
      END DO

!          BE10          y(j,75)
      DO j = 1, asize
        p(j) = em(j,I_BE10)
        l(j) = 0. + dd(j,I_BE10) + dw(j,I_BE10)
        y(j,I_BE10) = (yp(j,I_BE10)+dts*p(j))/(1.0+dts*l(j))
      END DO

!          RN222         y(j,76)
      DO j = 1, asize
        p(j) = em(j,I_RN222)
        l(j) = 0. + dd(j,I_RN222) + dw(j,I_RN222)                       &
     & +rc(j,291)
        y(j,I_RN222) = (yp(j,I_RN222)+dts*p(j))/(1.0+dts*l(j))
      END DO

!          PB210         y(j,77)
      DO j = 1, asize
        p(j) = em(j,I_PB210)                                            &
     & +rc(j,291)*y(j,I_RN222)
        l(j) = 0. + dd(j,I_PB210) + dw(j,I_PB210)
        y(j,I_PB210) = (yp(j,I_PB210)+dts*p(j))/(1.0+dts*l(j))
      END DO

!          STRAT         y(j,78)
      DO j = 1, asize
        p(j) = em(j,I_STRAT)
        l(j) = dd(j,I_STRAT)
        y(j,I_STRAT) = (yp(j,I_STRAT)+dts*p(j))/(1.0+dts*l(j))
      END DO

!         TROP          y(j,79)
      DO j = 1, asize
        p(j) = em(j,I_TROP)
        l(j) = dd(j,I_TROP)
        y(j,I_TROP) = (yp(j,I_TROP)+dts*p(j))/(1.0+dts*l(j))
      END DO

!       LNOX          Y(:,??)
! Lightning NOx produced could be in any NO-containg compound as cycling
! between NOz species is rapid. Losses below are permanant NOz sinks, sc
! as [LNOx]/Total [NOz]. HO2NO2 not included as has no permanant sinks.
!       P = EM(??)
!       L = (DD(5) + DW(5) + DD(13) + DW(13) + DD(58) + DW(58)
!    &    + DD(69) + DW(69)) * Y(??) /
!    &    (Y(4)+Y(5)+Y(6)+Y(7)+Y(13)+Y(55)+Y(56)+Y(57)Y(58)+Y(69))

          END DO  ! End of iteration loop
!
! At end of iteration, calculate flux terms.
! Loop over chunk of air parcels. Need to calculate fluxes for each
! air parcel for the cellflux calculation.
! For vectorization the following loop has to be split.
! First dimension chunk in the flux array keeps the data
! for all points in a chunk and allows splitting of the loop
! at arbitrary locations

! Set flux array to zero
          flux = 0.0

          DO i = 1, asize
!
            dtsm = dts / m(i)

! Calculate all fluxes
!
!    -----------------
!    THERMAL REACTIONS
!    -----------------
!
!      O3P + O2 + M = O3 + M
            flux(i,1)=flux(i,1)+rc(i,1)*y(i,2)*dtsm
!      O3P + NO + M = NO2 + M
            flux(i,5)=flux(i,5)+rc(i,5)*y(i,2)*y(i,4)*dtsm
!      O1D + M = O3P + M
            flux(i,7)=flux(i,7)+rc(i,7)*y(i,1)*dtsm
!      O1D + H2O = OH + OH
            flux(i,8)=flux(i,8)+rc(i,8)*y(i,1)*y(i,73)*dtsm
!      NO + O3 = NO2 + O2
            flux(i,11)=flux(i,11)+rc(i,11)*y(i,4)*y(i,11)*dtsm
!      NO2 + O3 = NO3 + O2
            flux(i,12)=flux(i,12)+rc(i,12)*y(i,5)*y(i,11)*dtsm
!      OH + O3 = HO2 + O2
            flux(i,13)=flux(i,13)+rc(i,13)*y(i,3)*y(i,11)*dtsm
!      HO2 + O3 = OH + 2O2
            flux(i,14)=flux(i,14)+rc(i,14)*y(i,16)*y(i,11)*dtsm
!      NO + NO3 = 2NO2
            flux(i,15)=flux(i,15)+rc(i,15)*y(i,4)*y(i,6)*dtsm
!      NO2 + O3P = NO
            flux(i,16)=flux(i,16)+rc(i,16)*y(i,2)*y(i,5)*dtsm
!      NO + HO2 = OH + NO2
            flux(i,17)=flux(i,17)+rc(i,17)*y(i,4)*y(i,16)*dtsm
!      NO2 + NO3 = NO + NO2 + O2
            flux(i,19)=flux(i,19)+rc(i,19)*y(i,5)*y(i,6)*dtsm
!      NO2 + NO3 = N2O5
            flux(i,20)=flux(i,20)+rc(i,20)*y(i,5)*y(i,6)*dtsm
!      NO2 + OH = HNO3
            flux(i,21)=flux(i,21)+rc(i,21)*y(i,5)*y(i,3)*dtsm
!      NO2 + HO2 = HO2NO2
            flux(i,22)=flux(i,22)+rc(i,22)*y(i,5)*y(i,16)*dtsm
!      HO2NO2 + M = NO2 + HO2
            flux(i,23)=flux(i,23)+rc(i,23)*y(i,65)*dtsm
!      HO2NO2 + OH = NO2
            flux(i,24)=flux(i,24)+rc(i,24)*y(i,65)*y(i,3)*dtsm
!      NO3 + NO3 = 2NO2
            flux(i,27)=flux(i,27)+rc(i,27)*y(i,6)*y(i,6)*dtsm
!      N2O5 = NO2 + NO3
            flux(i,29)=flux(i,29)+rc(i,29)*y(i,7)*dtsm
!      HO2 + OH = H2O + O2
            flux(i,30)=flux(i,30)+rc(i,30)*y(i,16)*y(i,3)*dtsm
!      OH + H2O2 = HO2 + H2O
            flux(i,31)=flux(i,31)+rc(i,31)*y(i,3)*y(i,14)*dtsm
!      NO3 + HO2 = HNO3
            flux(i,32)=flux(i,32)+rc(i,32)*y(i,6)*y(i,16)*dtsm
!      OH + H2 = HO2 + H2O
            flux(i,33)=flux(i,33)+rc(i,33)*y(i,3)*y(i,12)*dtsm
!      NO3 + HO2 = NO2 + OH
            flux(i,34)=flux(i,34)+rc(i,34)*y(i,6)*y(i,16)*dtsm
!      OH + HNO3 = NO3 + H2O
            flux(i,35)=flux(i,35)+rc(i,35)*y(i,3)*y(i,13)*dtsm
!      HO2 + HO2 = H2O2 + O2
            flux(i,36)=flux(i,36)+rc(i,36)*y(i,16)*y(i,16)*dtsm
!      OH + SO2 = HO2 + SA
            flux(i,39)=flux(i,39)+rc(i,39)*y(i,3)*y(i,26)*dtsm
!      N2O5 = 2NAER
            flux(i,41)=flux(i,41)+rc(i,43)*y(i,7)*dtsm
!      HNO3 = NAER
            flux(i,42)=flux(i,42)+rc(i,43)*y(i,13)*dtsm
!      ORGNIT = NAER
            flux(i,43)=flux(i,43)+rc(i,43)*y(i,69)*dtsm
!      CH3OOH + OH = CH3O2
            flux(i,44)=flux(i,44)+rc(i,44)*y(i,22)*y(i,3)*dtsm
!      CH3OOH + OH = HCHO + OH
            flux(i,45)=flux(i,45)+rc(i,45)*y(i,22)*y(i,3)*dtsm
!      C2H5OOH + OH = C2H5O2
            flux(i,46)=flux(i,46)+rc(i,44)*y(i,34)*y(i,3)*dtsm
!      C3H7OOH + OH = C3H7O2
            flux(i,47)=flux(i,47)+rc(i,44)*y(i,33)*y(i,3)*dtsm
!      C4H9OOH + OH = SC4H9O2
            flux(i,48)=flux(i,48)+rc(i,44)*y(i,35)*y(i,3)*dtsm
!      CH4 + OH = CH3O2 + H2O
            flux(i,59)=flux(i,59)+rc(i,59)*y(i,9)*y(i,3)*dtsm
!      NO + CH3O2 = HCHO + HO2 + NO2
            flux(i,60)=flux(i,60)+rc(i,60)*y(i,4)*y(i,15)*dtsm
!      CH3O2 + CH3O2 = 2HCHO + 2HO2
            flux(i,61)=flux(i,61)+rc(i,61)*y(i,15)*y(i,15)*dtsm
!      CH3O2 + CH3O2 = 2HCHO
            flux(i,62)=flux(i,62)+rc(i,62)*y(i,15)*y(i,15)*dtsm
!      CH3OH + OH = HCHO + HO2
            flux(i,63)=flux(i,63)+rc(i,63)*y(i,36)*y(i,3)*dtsm
!      CH3O2 + HO2 = CH3OOH + O2
            flux(i,65)=flux(i,65)+rc(i,65)*y(i,15)*y(i,16)*dtsm
!      HCHO + OH = HO2 + CO + H2O
            flux(i,66)=flux(i,66)+rc(i,66)*y(i,10)*y(i,3)*dtsm
!      HCHO + NO3 = HO2 + CO + HNO3
            flux(i,67)=flux(i,67)+rc(i,67)*y(i,10)*y(i,6)*dtsm
!      CO + OH = HO2 + CO2
            flux(i,70)=flux(i,70)+rc(i,70)*y(i,8)*y(i,3)*dtsm
!      OH + C2H6 = C2H5O2
            flux(i,71)=flux(i,71)+rc(i,71)*y(i,3)*y(i,17)*dtsm
!      C2H5O2 + NO = CH3CHO + NO2 + HO2
            flux(i,72)=flux(i,72)+rc(i,72)*y(i,18)*y(i,4)*dtsm
!      C2H5O2 + CH3O2 = CH3CHO + HCHO + HO2
            flux(i,73)=flux(i,73)+rc(i,73)*y(i,18)*y(i,15)*dtsm
!      CH3COO2 + CH3O2 = HCHO
            flux(i,74)=flux(i,74)+rc(i,74)*y(i,20)*y(i,15)*dtsm
!      OH + CH3CHO = CH3COO2
            flux(i,75)=flux(i,75)+rc(i,75)*y(i,3)*y(i,19)*dtsm
!      CH3COO2 + NO2 = PAN
            flux(i,77)=flux(i,77)+rc(i,77)*y(i,20)*y(i,5)*dtsm
!       PAN     +M       =CH3COO2 +NO2
            flux(i,78)=flux(i,78)+rc(i,78)*y(i,21)*dtsm
!       CH3COO2 +NO      =CH3O2   +NO2
            flux(i,79)=flux(i,79)+rc(i,79)*y(i,20)*y(i,4)*dtsm
!       CH3O2   +CH3COO2 =HCHO    +HO2     +CH3O2
            flux(i,80)=flux(i,80)+rc(i,80)*y(i,15)*y(i,20)*dtsm
!       OH      +NC4H10  =SC4H9O2
            flux(i,81)=flux(i,81)+rc(i,81)*y(i,3)*y(i,23)*dtsm
!       SC4H9O2 +NO      =CH3COE  +NO2     +HO2
            flux(i,83)=flux(i,83)+rc(i,83)*y(i,24)*y(i,4)*dtsm
!       SC4H9O2 +CH3O2   =CH3COE  +HO2     +HCHO
            flux(i,84)=flux(i,84)+rc(i,84)*y(i,24)*y(i,15)*dtsm
!       OH      +CH3COE  =CH3COX
            flux(i,86)=flux(i,86)+rc(i,86)*y(i,3)*y(i,25)*dtsm
!       C2H5O2  +C2H5O2  =CH3CHO  +HO2
            flux(i,90)=flux(i,90)+rc(i,90)*y(i,18)*y(i,18)*dtsm
!       CH3COO2 + CH3COO2 = 2.0*CH3O2
            flux(i,91)=flux(i,91)+rc(i,91)*y(i,20)*y(i,20)*dtsm
!       C3H8 + OH = C3H7O2
            flux(i,92)=flux(i,92)+rc(i,92)*y(i,31)*y(i,3)*dtsm
!       C3H7O2 + NO = NO2 + HO2 + ACETONE
            flux(i,93)=flux(i,93)+rc(i,93)*y(i,32)*y(i,4)*dtsm
!       ACETONE + OH = ACETO2
            flux(i,94)=flux(i,94)+rc(i,94)*y(i,37)*y(i,3)*dtsm
!       ACETO2 + NO = NO2 + HCHO + CH3COO2
            flux(i,95)=flux(i,95)+rc(i,95)*y(i,38)*y(i,4)*dtsm
!       ACETO2 + CH3O2 = HO2 + HCHO + CH3COO2
            flux(i,96)=flux(i,96)+rc(i,96)*y(i,38)*y(i,15)*dtsm
!       CH37O2 + CH3O2 = HO2 + HCHO + ACETONE
            flux(i,97)=flux(i,97)+rc(i,97)*y(i,32)*y(i,15)*dtsm
!       PAN + OH = NO3 + HCHO
            flux(i,98)=flux(i,98)+rc(i,98)*y(i,21)*y(i,3)*dtsm
!       C2H5O2 + HO2 = C2H5OOH
            flux(i,99)=flux(i,99)+rc(i,99)*y(i,18)*y(i,16)*dtsm
!       C2H5OOH + OH = CH3CHO + OH
            flux(i,100)=flux(i,100)+rc(i,100)*y(i,34)*y(i,3)*dtsm
!       C3H7O2 + HO2 = C3H7OOH
            flux(i,101)=flux(i,101)+rc(i,101)*y(i,32)*y(i,16)*dtsm
!       C3H7OOH + OH = ACETONE + OH
            flux(i,102)=flux(i,102)+rc(i,102)*y(i,33)*y(i,3)*dtsm
!       SC4H9O2 + HO2 = C4H9OOH
            flux(i,103)=flux(i,103)+rc(i,103)*y(i,24)*y(i,16)*dtsm
!       C4H9OOH + OH = CH3COE + OH
            flux(i,104)=flux(i,104)+rc(i,104)*y(i,35)*y(i,3)*dtsm
!        NO + CH3COX = CH3COO2 + NO2 + CH3CHO
            flux(i,105)=flux(i,105)+rc(i,105)*y(i,40)*y(i,4)*dtsm
!        CH3O2 + CH3COX = HCHO + 2.0*HO2 + CH3COY
            flux(i,106)=flux(i,106)+rc(i,106)*y(i,15)*y(i,40)*dtsm
!        OH      +C2H4    =CH2O2C
            flux(i,109)=flux(i,109)+rc(i,109)*y(i,3)*y(i,27)*dtsm
!        CH2O2C +NO = 2.0*HCHO + HO2 + NO2
            flux(i,110)=flux(i,110)+rc(i,110)*y(i,41)*y(i,4)*dtsm
!        CH3O2   +CH2O2C  = 3.0*HCHO + 2.0*HO2
            flux(i,111)=flux(i,111)+rc(i,111)*y(i,15)*y(i,41)*dtsm
!        O3 + C2H4 = HCHO + 0.31*CO + 0.13*H2 +0.20*HO2 +0.47CH2OO
            flux(i,112)=flux(i,112)+rc(i,112)*y(i,11)*y(i,27)*dtsm
!        O3 + C3H6 = HCHO +.3*CH4 +.4*CO +.28*OH +.30*HO2 +.58*CH3O2 +.1
            flux(i,123)=flux(i,123)+rc(i,123)*y(i,11)*y(i,28)*dtsm
!        O3 + C3H6 = CH3CHO + 0.24*H2 + 0.58*CO +0.18HO2
            flux(i,124)=flux(i,124)+rc(i,124)*y(i,11)*y(i,28)*dtsm
!        OH + C3H6 = CH3CHX
            flux(i,125)=flux(i,125)+rc(i,125)*y(i,3)*y(i,28)*dtsm
!        CH3CHX + NO = HCHO + HO2 + CH3CHO + NO2
            flux(i,126)=flux(i,126)+rc(i,126)*y(i,43)*y(i,4)*dtsm
!        CH3O2 + CH3CHX = 2.0*HCHO + 2.0*HO2 + CH3CHO
            flux(i,127)=flux(i,127)+rc(i,127)*y(i,15)*y(i,43)*dtsm
!        C5H8 + O3 = MVK + .78CO + .22CH2OO + .27HO2 + .27OH
            flux(i,128)=flux(i,128)+rc(i,128)*y(i,48)*y(i,11)*dtsm
!        MVK + O3 = MGLYOX + .76CO + .24CH2OO + .36HO2 + .36OH
            flux(i,129)=flux(i,129)+rc(i,129)*y(i,50)*y(i,11)*dtsm
!        CH2OO + NO = NO2 + HCHO
!           flux(i,130)=flux(i,130)+rc(i,130)*y(i,39)*y(i,4)*dtsm
!        CH2OO + NO2 = NO3 + HCHO
!           flux(i,131)=flux(i,131)+rc(i,131)*y(i,39)*y(i,5)*dtsm
!        CH2OO =  HCHO
!           flux(i,132)=flux(i,132)+rc(i,132)*y(i,39)*dtsm
!        CH2OO + SO2 =  SA + HCHO
!           flux(i,133)=flux(i,133)+rc(i,133)*y(i,39)*y(i,26)*dtsm
!        NO3 + C2H6 = C2H5O2 + HNO3
            flux(i,201)=flux(i,201)+rc(i,201)*y(i,6)*y(i,17)*dtsm
!        NO3 + NC4H10 = SC4H9O2 + HNO3
            flux(i,202)=flux(i,202)+rc(i,202)*y(i,6)*y(i,23)*dtsm
!        NO3 + C2H4 = RNC2H4
            flux(i,203)=flux(i,203)+rc(i,203)*y(i,6)*y(i,27)*dtsm
!        RNC2H4 + OH = HCHO + NO2
            flux(i,204)=flux(i,204)+rc(i,204)*y(i,55)*y(i,3)*dtsm
!        NO3 + C3H6 = RNC3H6
            flux(i,205)=flux(i,205)+rc(i,205)*y(i,6)*y(i,28)*dtsm
!        RNC3H6 + OH = CH3CHO + NO2
            flux(i,206)=flux(i,206)+rc(i,206)*y(i,56)*y(i,3)*dtsm
!        NO3 + CH3CHO = CH3COO2 + HNO3
            flux(i,208)=flux(i,208)+rc(i,208)*y(i,6)*y(i,19)*dtsm
!        NO3 + C5H8 = RNC5H8 + HO2
            flux(i,209)=flux(i,209)+rc(i,209)*y(i,6)*y(i,48)*dtsm
!        RNC5H8 + OH = MVK + NO2
            flux(i,210)=flux(i,210)+rc(i,210)*y(i,57)*y(i,3)*dtsm
!        HO2 + TOLP1 = HCHO + HO2 + MEMALD + GLYOX
            flux(i,212)=flux(i,212)+rc(i,212)*y(i,16)*y(i,59)*dtsm
!        NO2 + TOLP1 = MEMALD + GLYOX + HO2 + NO2
            flux(i,213)=flux(i,213)+rc(i,213)*y(i,5)*y(i,59)*dtsm
!        OH + DMS = CH3SO + HCHO
            flux(i,214)=flux(i,214)+rc(i,214)*y(i,3)*y(i,60)*dtsm
!        OH + DMS = DMSO + HO2
            flux(i,215)=flux(i,215)+rc(i,215)*y(i,3)*y(i,60)*dtsm
!        CH3SO + O3 = CH3SO2
            flux(i,216)=flux(i,216)+rc(i,216)*y(i,66)*y(i,11)*dtsm
!        CH3SO + NO2 = CH3SO2 + NO
            flux(i,217)=flux(i,217)+rc(i,217)*y(i,66)*y(i,5)*dtsm
!        CH3SO2 + O3 = CH3SO3
            flux(i,218)=flux(i,218)+rc(i,218)*y(i,64)*y(i,11)*dtsm
!        CH3SO2 + NO2 = CH3SO3 + NO
            flux(i,219)=flux(i,219)+rc(i,219)*y(i,64)*y(i,5)*dtsm
!        CH3SO2 + O2 = CH3O2 + SO2
            flux(i,220)=flux(i,220)+rc(i,220)*y(i,64)*dtsm
!        CH3SO3 + HO2 = MSA
            flux(i,221)=flux(i,221)+rc(i,221)*y(i,67)*y(i,16)*dtsm
!        CH3SO3 + O2 = CH3O2 + SA
            flux(i,222)=flux(i,222)+rc(i,222)*y(i,67)*dtsm
!        CH3SO3 + HCHO = MSA + HO2 + CO
            flux(i,223)=flux(i,223)+rc(i,223)*y(i,67)*y(i,10)*dtsm
!        OH + DMSO = DMSO2 + HO2
            flux(i,224)=flux(i,224)+rc(i,224)*y(i,3)*y(i,62)*dtsm
!        OH + DMSO2 = DMSP
            flux(i,225)=flux(i,225)+rc(i,225)*y(i,3)*y(i,61)*dtsm
!        DMSP + NO = NO2 + HCHO + CH3SO2
            flux(i,226)=flux(i,226)+rc(i,226)*y(i,4)*y(i,63)*dtsm
!        DMSP + CH3O2 = HO2 + 2HCHO + CH3SO2
            flux(i,227)=flux(i,227)+rc(i,227)*y(i,63)*y(i,15)*dtsm
!        DMS + NO3 = CH3SO + HCHO + HNO3
            flux(i,228)=flux(i,228)+rc(i,228)*y(i,60)*y(i,6)*dtsm
!        OH + OXYL = HO2 + MEMALD + MGLYOX
            flux(i,230)=flux(i,230)+rc(i,230)*y(i,3)*y(i,29)*dtsm
!        OXYL1 + NO2 = MEMALD + NO2 + HO2 + MGLYOX
            flux(i,231)=flux(i,231)+rc(i,231)*y(i,45)*y(i,5)*dtsm
!        OH + MEMALD = MEMALD1
            flux(i,232)=flux(i,232)+rc(i,232)*y(i,3)*y(i,46)*dtsm
!        MEMALD1 + NO = HO2 + NO2 + GLYOX + MGLYOX
            flux(i,233)=flux(i,233)+rc(i,233)*y(i,47)*y(i,4)*dtsm
!        OH + TOLUEN = HO2 + MEMALD + GLYOX
            flux(i,234)=flux(i,234)+rc(i,234)*y(i,3)*y(i,54)*dtsm
!        OH + TOLUEN = TOLP1
            flux(i,235)=flux(i,235)+rc(i,235)*y(i,3)*y(i,54)*dtsm
!        OH + OXYL = OXYL1
            flux(i,236)=flux(i,236)+rc(i,236)*y(i,3)*y(i,29)*dtsm
!        OH + ORGNIT = MEMALD + GLYOX + NO2
            flux(i,237)=flux(i,237)+rc(i,237)*y(i,3)*y(i,69)*dtsm
!        NO3 + ORGNIT = MEMALD + GLYOX + NO2
            flux(i,238)=flux(i,238)+rc(i,238)*y(i,6)*y(i,69)*dtsm
!        OXYL1 + HO2 = MGLYOX + MEMALD
            flux(i,240)=flux(i,240)+rc(i,240)*y(i,45)*y(i,16)*dtsm
!        MEMALD1 + CH3O2 = 2.0*HO2 + HCHO + GLYOX + MGLYOX
            flux(i,241)=flux(i,241)+rc(i,241)*y(i,47)*y(i,15)*dtsm
!        RO2IP1 + CH3O2 = 2.0*HO2 + 2.0*HCHO + MVK
!           flux(i,242)=flux(i,242)+rc(i,242)*y(i,49)*y(i,15)*dtsm
!        RO2IP2 + CH3O2 = 2.0*HO2 + 2.0*HCHO + MGLYOX
!           flux(i,243)=flux(i,243)+rc(i,243)*y(i,51)*y(i,15)*dtsm
!        RO2IP1 + HO2 = ISOPOOH
            flux(i,244)=flux(i,244)+rc(i,244)*y(i,49)*y(i,16)*dtsm
!        ISOPOOH + OH = MVK + HCHO + OH
            flux(i,245)=flux(i,245)+rc(i,245)*y(i,52)*y(i,3)*dtsm
!        RO2IP2 + HO2 = MVKOOH
            flux(i,246)=flux(i,246)+rc(i,246)*y(i,51)*y(i,16)*dtsm
!        MVKOOH + OH = MGLYOX + HCHO + OH
            flux(i,247)=flux(i,247)+rc(i,247)*y(i,53)*y(i,3)*dtsm
!        MGLYOX + OH = CH3COO2 + CO
            flux(i,248)=flux(i,248)+rc(i,248)*y(i,42)*y(i,3)*dtsm
!        GLYOX + OH = HO2 + CO
            flux(i,249)=flux(i,249)+rc(i,249)*y(i,44)*y(i,3)*dtsm
!        OH + C5H8 = RO2IP1
            flux(i,251)=flux(i,251)+rc(i,251)*y(i,3)*y(i,48)*dtsm
!        RO2IP1 + NO = MVK + HO2 + HCHO + NO2
            flux(i,252)=flux(i,252)+rc(i,252)*y(i,49)*y(i,4)*dtsm
!        OH + MVK = RO2IP2
            flux(i,253)=flux(i,253)+rc(i,253)*y(i,3)*y(i,50)*dtsm
!        RO2IP2 + NO = MGLYOX + HCHO + HO2 + NO2
            flux(i,254)=flux(i,254)+rc(i,254)*y(i,51)*y(i,4)*dtsm
!        HSO3-(aq) + H2O2(aq) = SO4--(aq)
            flux(i,260)=flux(i,260)+faq(i,2)*dtsm
!        HSO3-(aq) + O3(aq) = SO4--(aq)
            flux(i,261)=flux(i,261)+faq(i,4)*dtsm
!        SO3--(aq) + O3(aq) = SO4--(aq)
            flux(i,262)=flux(i,262)+faq(i,5)*dtsm
!        SO2(aq) = SO4--(aq)
            flux(i,263)=flux(i,263)+faq(i,1)*dtsm
!        2NH4+(aq) + SO4--(aq) = (NH4)2SO4 (aq)
            flux(i,264)=flux(i,264)+faq(i,3)*dtsm

!      ------------------
!      RADIOACTIVE DECAYS
!      ------------------
!
!        Be-7 =
            flux(i,290)=flux(i,290)+rc(i,290)*y(i,i_be7)*dtsm
!        Rn-222 = Pb-210
            flux(i,291)=flux(i,291)+rc(i,291)*y(i,i_rn222)*dtsm
!
!
!    --------------------
!    PHOTOLYTIC REACTIONS
!    --------------------
!
! O3 + hv = OP
            flux(i,301) = kj(i,1)*y(i,11)*dtsm
! O3 + hv = OD
            flux(i,302) = kj(i,2)*y(i,11)*dtsm
! NO2 + hv = NO + OP
            flux(i,303) = kj(i,3)*y(i,5)*dtsm
! H2O2 + hv = OH + OH
            flux(i,304) = kj(i,4)*y(i,14)*dtsm
! HNO3 + hv = NO2 + OH
            flux(i,305) = kj(i,5)*y(i,13)*dtsm
! HCHO + hv = CO + HO2 + HO2
            flux(i,306) = kj(i,6)*y(i,10)*dtsm
! HCHO + hv = CO + H2
            flux(i,307) = kj(i,7)*y(i,10)*dtsm
! CH3CHO + hv = CH3O2 + HO2 + CO
            flux(i,308) = kj(i,8)*y(i,19)*dtsm
! CH3COE + hv = C2H5O2 + CH3COO2
            flux(i,309) = kj(i,9)*y(i,25)*dtsm
! HO2NO2 + hv = HO2 + NO2
            flux(i,310) = kj(i,10)*y(i,65)*dtsm
! MGLYOX + hv = CH3COO2 + HO2 + CO
            flux(i,311) = kj(i,11)*y(i,42)*dtsm
! GLYOX + hv  = HCHO + CO
            flux(i,312) = kj(i,12)*y(i,44)*dtsm
! NO3 + hv = NO
            flux(i,313) = kj(i,13)*y(i,6)*dtsm
! NO3 + hv = NO2 + OP
            flux(i,314) = kj(i,14)*y(i,6)*dtsm
! N2O5 + hv = NO2 + NO3
            flux(i,315) = kj(i,15)*y(i,7)*dtsm
! CH3OOH + hv = HCHO + HO2 + OH
            flux(i,316) = kj(i,16)*y(i,22)*dtsm
! PAN + hv = CH3COO2 + NO2
            flux(i,317) = kj(i,17)*y(i,21)*dtsm
! ACETONE + hv = CH3COO2 + CH3O2
            flux(i,318) = kj(i,9)*y(i,37)*dtsm
! C2H5OOH + hv = CH3CHO + HO2 + OH
            flux(i,319) = kj(i,16)*y(i,34)*dtsm
! C3H7OOH + hv = ACETONE + HO2 + OH
            flux(i,320) = kj(i,16)*y(i,33)*dtsm
! C4H9OOH + hv = CH3COE + HO2 + OH
            flux(i,321) = kj(i,16)*y(i,35)*dtsm
! ISOPOOH + hv = MVK + HO2 + OH + HCHO
            flux(i,322) = kj(i,16)*y(i,52)*dtsm
! MVKOOH + hv = MGLYOX + HO2 + OH + HCHO
            flux(i,323) = kj(i,16)*y(i,53)*dtsm
!
          END DO

!    -------------------------
!    EMISSIONS AND DEPOSITION:
!    -------------------------
!
          dtsm_a(1:asize) = dts / m(1:asize)
          DO kk = 1, nc
            flux(1:asize,400+kk) = em(1:asize,kk) * dtsm_a(1:asize)
            flux(1:asize,500+kk) = y(1:asize,kk) * dd(1:asize,kk)       &
     &        * dtsm_a(1:asize)
            flux(1:asize,600+kk) = y(1:asize,kk) * dw(1:asize,kk)       &
     &        * dtsm_a(1:asize)
          END DO

! Which value of m to use ??
          DO i = 1, asize
            flux(i,400+nc+1) = em(i,nc+1) * dts / m(1)
            flux(i,400+nc+2) = em(i,nc+2) * dts / m(1)
!
!  Calculate ozone budget terms:
!
!  HO2+NO
            flux(i,350)=flux(i,17)
!  CH3O2+NO
            flux(i,351)=flux(i,60)
!  RO2+NO
            flux(i,352)=                                                &
     &         +flux(i,72)                                              &
                                              ! C2H5O2+NO
     &         +flux(i,79)                                              &
                                              ! CH3COO2+NO
     &         +flux(i,83)                                              &
                                              ! SC4H9O2+NO
     &         +flux(i,93)                                              &
                                              ! C3H7O2+NO
     &         +flux(i,95)                                              &
                                              ! ACETO2+NO
     &         +flux(i,105)                                             &
                                              ! CH3COX+NO
     &         +flux(i,110)                                             &
                                              ! CH2O2C+NO
     &         +flux(i,126)                                             &
                                              ! CH3CHX+NO
     &         +flux(i,252)                                             &
                                              ! RO2IP1+NO
     &         +flux(i,254)                                             &
                                              ! RO2IP2+NO
     &         +flux(i,226)                                             &
                                              ! DMSP+NO
     &         +flux(i,233)                   ! MEMALD1+NO

!  Stratospheric input
            flux(i,353)=flux(i,411)

!  Total Production = (HO2+NO)+(CH3O2+NO)+(RO2+NO)+(Strat O3)
            flux(i,354)=flux(i,17)+flux(i,60)+flux(i,352)+flux(i,411)

!  Ozone destruction:
!  O1D+H2O
            flux(i,355)=flux(i,8)
!  O3+OH
            flux(i,356)=flux(i,13)
!  O3+HO2
            flux(i,357)=flux(i,14)
!  O3 + Hydrocarbons
            flux(i,358)=                                                &
     &         +flux(i,112)                                             &
                                              ! O3+C2H4 =>
     &         +flux(i,123)                                             &
                                              ! O3+C3H6 =>
     &         +flux(i,124)                                             &
                                              ! O3+C3H6 =>
     &         +flux(i,128)                                             &
                                              ! O3+C5H8 =>
     &         +flux(i,129)                                             &
                                              ! O3+MVK =>
     &         +flux(i,216)                                             &
                                              ! O3+CH3SO =>
     &         +flux(i,218)                   ! O3+CH3SO2 =>

!  Other net gains and losses
            flux(i,359)=                                                &
     &         -flux(i,305)                                             &
                                              ! HNO3+hv => OH+NO2
     &         -flux(i,314)*2.0                                         &
                                              ! NO3+hv => NO2+O3P
     &         -flux(i,315)                                             &
                                              ! N2O5+hv => NO2+NO3
     &         +flux(i,12)*2.0                                          &
                                              ! NO2+O3 => NO3+O2
     &         -flux(i,15)*2.0                                          &
                                              ! NO+NO3 => 2NO2
     &         +flux(i,16)*2.0                                          &
                                              ! NO2+O3P => NO+O2
     &         +flux(i,20)                                              &
                                              ! NO2+NO3 => N2O5
     &         +flux(i,21)                                              &
                                              ! NO2+OH => HNO3
     &         -flux(i,27)*2.0                                          &
                                              ! NO3+NO3 => 2NO2
     &         -flux(i,29)                                              &
                                              ! N2O5 => NO2+NO3
     &         -flux(i,34)                                              &
                                              ! NO3+HO2 => OH+NO2
     &         +flux(i,98)                                              &
                                              ! OH+PAN => NO3+HCHO
     &         -flux(i,204)                                             &
                                              ! OH+RNC2H4 => HCHO+NO2
     &         -flux(i,206)                                             &
                                              ! OH+RNC3H6 => CH3CHO+NO2
     &         -flux(i,210)                                             &
                                              ! OH+RNC5H8 => CH3CHO+NO2
     &         +flux(i,213)                                             &
                                              ! NO2+TOLP1 => ORGNIT
     &         +flux(i,217)                                             &
                                              ! NO2+CH3SO => CH3SO2+NO
     &         +flux(i,219)                                             &
                                              ! NO2+CH3SO2 => CH3SO3+NO
     &         +flux(i,231)                                             &
                                              ! NO2+OXYL1 => ORGNIT
     &         +flux(i,505)                                             &
                                              ! NO2 dry deposition
     &         +flux(i,521)                   ! PAN dry deposition
!  O3 Dry Deposition
            flux(i,360)=flux(i,511)
!  Total Destruction = (O1D+H2O)+(O3+OH)+(O3+HO2)+
!                      (other)+(O3+HCs)+(O3 dry depn.)
            flux(i,361)=                                                &
     &         +flux(i,8)+flux(i,13)+flux(i,14)                         &
     &         +flux(i,359)+flux(i,358)+flux(i,511)
!  Net Production
            flux(i,362)=flux(i,354)-flux(i,361)
!  Net Chemical Production
            flux(i,363)=flux(i,362)-flux(i,411)+flux(i,511)
!
          END DO ! DO i = 1, asize
!
! Store total fluxes and 3d fluxes
          DO kk = 1, nflux
            k1 = flist(1,kk)
            k2 = flist(2,kk)
            IF (k1 == 1799) THEN
!
! Calculate AOT40. Calculate it everywhere for now. Units at this
! stage are ppbv-s.
              DO i = 1, asize
                flux(i,700) = MAX((1.0e9*y(i,i_o3) / m(i))-40.0,0.0)
                totflu(kk) = totflu(kk) + flux(i,700)
              END DO
              IF (k2 > 0) THEN
                DO i = 1, asize
                  j = i + j0 - 1
                  cellflux(k2,j) = cellflux(k2,j) + flux(i,700)
                END DO
              END IF
            END IF
          END DO
!
! Store total fluxes and 3d fluxes
          DO kk = 1, nflux
            k1 = flist(1,kk)
            k2 = flist(2,kk)
            IF (k1 < 700) THEN
              IF (k2 > 0) THEN
                DO i = 1, asize
                  j = i + j0 - 1
                  totflu(kk) = totflu(kk) + flux(i,k1)
                  cellflux(k2,j) = cellflux(k2,j) + flux(i,k1)
                END DO
              ELSE
                DO i = 1, asize
                  totflu(kk) = totflu(kk) + flux(i,k1)
                END DO
              END IF
            END IF
          END DO
!
          yp = y
!
! flux(i,1-NR) holds fluxes through reactions 1-NR
! flux(i,301-300+NDJ) holds fluxes through photolytic reactions 1-NDJ
! flux(i,350-356)    holds ozone budget fluxes
! flux(i,401-400+NC) holds fluxes through emissions of species 1-NC
! flux(i,501-500+NC) holds fluxes through depositions of species 1-NC
! flux(i,601-600+NC) holds fluxes through wet deposition of species 1-NC
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        END DO                    ! End of Lagrangian time integration

! Mix in the changes to H2O2,SO2,SA,NH4 and NH3.
        DO j = 1, asize
          IF (ll(j) > l_crit) THEN

! in gas phase  Y(14)=Y(14)-(C1H2O2-CTH2O2)*CVF*mol2molec
! in gas phase  Y(30)=Y(30)-(C1SO4-CS(10))*CVF*NA/1.0E+3

            y(j,26) = y(j,26) - (c1so2(j)-ctso2(j)) * mol2molec * cvf
            y(j,39) = y(j,39) - (c1nh4(j)-ctnh4(j)) * mol2molec * cvf
            y(j,70) = y(j,70) - (c1nh3(j)-ctnh3(j)) * mol2molec * cvf
            IF (y(j,26) < 0.0) y(j,26) = m(j) * small_conc
            IF (y(j,39) < 0.0) y(j,39) = m(j) * small_conc
            IF (y(j,70) < 0.0) y(j,70) = m(j) * small_conc
          END IF
        END DO

! Convert final concentrations back to mixing ratio
        DO i = 1, nc
          xx(i,j0:j1) = y(1:asize,i) / m(1:asize)
        END DO

      END DO  ! End of loop over cells  do j0 = 1, nfill, chunk

      END SUBROUTINE DERIV
#endif
