#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE STOCHEM(                                               &
     &  u,v,q,sa,t0,bl,sice,orog,w,dc,p,lnp,o3um,                       &
     &  clw,acr,t,cc,                                                   &
     &  acp,adp,tropz,up,upd,cc_base_level,cc_top_level,                &
     &  land_fraction,                                                  &
     &  totpar, dirpar, t0tile, gsf, lai_ft, soilmc,                    &
     &  canwc, stcon, ra, rq, so4_vd,                                   &
     &  year,month,day,                                                 &
     &  StochDumpRequest,                                               &
     &  UMStepNo,fixhd12,stochem_call_step,z_top_of_model,              &
     &  first_constant_r_rho_level,o3mmr,ch4mmr)
! ----------------------------------------------------------------------
! Purpose:
! Calculate concentrations of tropospheric trace gases.
! Main routine.
!
! Method:
!
! TEST TRIVIAL DELTA
! Original Programmers: Bill, Collins, Dick Derwent, Colin Johnson
!                      and David Stevenson 1993-2000
!
! Current code owner: Colin Johnson
!
! History:
! Version   Date                    Comment
!  5.2    02/08/01  Created from vn4.5 code. C.E. Johnson
!  5.5    11/02/04  Now uses NEC random number generating routine
!                   FRANV. seed arrays changed. M.G. Sanderson
!  5.5    11/03/04  Vectorisation updates. K. Ketelsen.
! 12/09/03    6.0         Correct call to EREPORT.  Remove repeated
!  6.0    12/09/03  Corrections. C.E. Johnson
!  6.1    20/08/04  STOCHEM methane/ozone coupling. C.E. Johnson
!                         declaration of POS1, POS2.          Paul Dando
!  6.2    21/10/05  Replace GSYNC with SSYNC. P.Selwood.
!  6.2    20/04/05  Faster load balancing routines. R. Johanni.
!  6.2    02/03/06  Extensive additions and changes to use new dry
!                   deposition and interactive isoprene emission
!                   schemes. M.G. Sanderson
!
!VVV  V5.2.1  STOCHEM  2/8/01 -
! ----------------------------------------------------------------------
!
!
!      EULERIAN GRID IN DEGREES : CONC(NC,NLONG,NLPE,NLEV),
!      GRID BEGINS AT 90 N, 0 E.
!
!      LAGRANGIAN CELLS : XX(NC,NFILL), CENTRES AT: POS(4,NFILL)
!                         IPOS(3,NFILL) HAS LONG, LAT, and ETA
!                         INDICES
!
!-----------------------------------------------------------------------

!   NLONG: No. of Eulerian longitude grids; East from Greenwich
!   NLAT:  No. of Eulerian latitude grids;  North to South
!   NLEV:  No. of Eulerian levels;
!   NCELL: No. of Lagrangian cells (nominal);
!   DLAT:  Eulerian latitudinal grid (degrees);
!   DLONG: Eulerian longitudinal grid  (degrees);
!   NMETLAT,NMETLONG,NMETLEV: Meteorological grid dimensions;
!   DLATM,DLONGM: Meteorological grid spacing.

!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_OUT
      USE IN_STOCHEM_INTF
      USE KK_GCOM_MOD
      IMPLICIT NONE
!-----------------------------------------------------------------------

      INTEGER, INTENT(IN)    :: umstepno
      INTEGER, INTENT(IN)    :: fixhd12       ! UM version No.
      INTEGER, INTENT(IN)    :: first_constant_r_rho_level
      REAL,    INTENT(IN)    :: z_top_of_model
      REAL,    INTENT(IN)    :: stochem_call_step
      INTEGER, INTENT(INOUT) :: year
      INTEGER, INTENT(INOUT) :: month
      REAL,    INTENT(INOUT) :: day

! Meteorological data from diagnostics:
! Those marked with a * are derived data.
! Those marked with a # are modified in READPP.
! Other diagnostics are used within READPP.

! U Component of Wind, V Component of Wind,
! T Temperature, Q Specific Humidity, CLW Cloud liquid water,
! CC Convective cloud fraction, DC dynamic cloud fraction,
! P_TH Pressure at theta levels, FLC fractional land cover,
! UP Convective mass flux, UPD Detrainment,
! W Component of wind,
! P0 Surface Pressure, T0 Surface Temperature,
! TROPZ Tropopause height,
! BL Boundary Layer Depth (#), SI Sea Ice fraction,
! SA Snow amount, OROG Orography, HF Heat Flux,
! ADR Dynamic rain, ADS Dynamic snow,
! ACR Convective rain, ADS Convective snow,
! ACP Convective Precipitation (*), ADP Dynamic Precipitation (*).

      INTEGER,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: cc_base_level
      INTEGER,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: cc_top_level

      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,nmetlev)   :: t
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,nmetlev)   :: q
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,nmetlev)   :: clw
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,nmetlev)   :: cc
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,nmetlev)   :: dc
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,0:nmetlev) :: p
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,0:nmetlev) :: lnp
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,nmetlev)   :: up
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,nmetlev)   :: upd
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,0:nmetlev) :: u
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,0:nmetlev) :: v
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,0:nmetlev) :: w
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,nmetlev)   :: o3um

      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe)           :: totpar
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe)           :: dirpar
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,ntype)     :: t0tile
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,ntype)     :: gsf
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,npft)      :: lai_ft
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,npft)      :: canwc
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe)           :: soilmc
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,npft)      :: stcon
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,ntype)     :: ra
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe,nc)        :: rq
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe)           :: so4_vd

      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: t0
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: bl
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: sice
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: sa
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: tropz
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: acr
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: land_fraction
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: acp
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: adp
      REAL,INTENT(IN),DIMENSION(nlonpe,nlatpe) :: orog
      REAL,INTENT(OUT),DIMENSION(nlonpe,nlatpe,nmetlev) :: o3mmr
      REAL,INTENT(OUT),DIMENSION(nlonpe,nlatpe,nmetlev) :: ch4mmr

      LOGICAL,       INTENT(INOUT)  :: StochDumpRequest
                                        ! likely
      INTEGER :: i       ! multi-use loop index
      INTEGER :: j       ! loop index over parcel numbers
      INTEGER :: k       ! loop index over vertical dimension
      INTEGER :: l       ! loop index over longitude
      INTEGER :: n       ! loop index over fluxes
      INTEGER :: k2      ! 3D flux number
      INTEGER, SAVE :: nchem    ! no. of 3D chemical fields requested
      INTEGER, SAVE :: nstation ! no. of stations for output
      INTEGER :: imsg=1  ! imsg is needed for some GCOM routines
      INTEGER :: info    ! info is needed for some GCOM routines
      INTEGER :: cellbase    ! used for indexing cell nos on each pe
      INTEGER, SAVE :: nfill ! no of cells currently carried by this pe
      INTEGER, SAVE :: nstat !
      INTEGER, SAVE :: nflux
      INTEGER, SAVE :: navg
      INTEGER, DIMENSION(:,:), ALLOCATABLE         :: cellswap
      INTEGER, DIMENSION(0:nproc-1)                :: max_cells
      INTEGER, DIMENSION(:,:), ALLOCATABLE         :: ipos
      INTEGER, DIMENSION(:,:), ALLOCATABLE         :: nbl
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: nm
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE       :: nnn
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: np
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE       :: nph
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: cellno
      INTEGER, DIMENSION(2,numflux),      SAVE :: flist
      INTEGER, DIMENSION(numchem),        SAVE :: clist
      INTEGER, DIMENSION(numflux),        SAVE :: stash_code
      INTEGER, DIMENSION(numfol),         SAVE :: follist
! each PE has own seeds
      INTEGER, DIMENSION(0:maxproc-1,ransize), SAVE :: seed2
      INTEGER, DIMENSION(0:maxproc-1,ransize), SAVE :: seed3

      INTEGER :: n3dflux
      INTEGER :: nmap
      INTEGER, DIMENSION(numflux) :: idxflux
      INTEGER, DIMENSION(0:nproc-1) :: nc_actual, nc_target
      INTEGER, DIMENSION(3,2*nproc) :: balance_map

      REAL, SAVE :: daym0
      REAL       :: time
      REAL, SAVE :: endtime
      REAL       :: period=0.0
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: pos
      REAL, DIMENSION(:,:), ALLOCATABLE       :: posold

!   Lagrangian mixing ratio: XX, Position: POS, Velocity: V,
!   Temperature: TL, & Random numbers: B.
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE       :: dja
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: dj
      REAL, DIMENSION(:,:,:),   ALLOCATABLE       :: cloud
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: o3_ls
      REAL, DIMENSION(:,:),     ALLOCATABLE, SAVE :: land
      REAL, DIMENSION(:,:),     ALLOCATABLE, SAVE :: o3col
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: seaicefr
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: snowfr
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: psurf
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: tsurf
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: tropdat
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: landum
      REAL, DIMENSION(:),       ALLOCATABLE       :: tl
      REAL, DIMENSION(:),       ALLOCATABLE       :: ql
      REAL, DIMENSION(:),       ALLOCATABLE       :: rhl
      REAL, DIMENSION(:),       ALLOCATABLE       :: m
      REAL, DIMENSION(:),       ALLOCATABLE       :: liq
      REAL, DIMENSION(:,:),     ALLOCATABLE, SAVE :: xx
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: em
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: dd
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: dw
      REAL, DIMENSION(:,:,:),   ALLOCATABLE       :: rh
      REAL, DIMENSION(:),       ALLOCATABLE       :: qs

      REAL, DIMENSION(:,:),     ALLOCATABLE       :: terpem

      REAL, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: m0
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE       :: mass
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE       :: m1
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE       :: rc
      REAL, DIMENSION(:,:),     ALLOCATABLE       :: o3_stom_frac

! Eulerian concentration: CONC,
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: conc
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: mconc
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: sdconc
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: o3conc
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: ch4conc
      REAL, DIMENSION(2,nc),                 SAVE :: totm0
      REAL, DIMENSION(2,nc),                 SAVE :: totavg
      REAL, DIMENSION(2,nc),                 SAVE :: totmas
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: tflux
      REAL, DIMENSION(:,:),     ALLOCATABLE, SAVE :: cellflux
      REAL, DIMENSION(numflux),              SAVE :: totflu
      REAL, DIMENSION(numstat),              SAVE :: stlon
      REAL, DIMENSION(numstat),              SAVE :: stlat
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: emiss
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: estore
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: dstore
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: so2em
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: be7em
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: be10em
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: lnoxem
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: acnoxem
      REAL, DIMENSION(:,:,:),   ALLOCATABLE, SAVE :: nflashes
      REAL,                                  SAVE :: totso2em
      REAL,                                  SAVE :: totbe7em
      REAL,                                  SAVE :: totbe10em
      REAL,                                  SAVE :: totlnoxem
      REAL,                                  SAVE :: totacnoxem
      REAL :: mtemp
      REAL, DIMENSION(nlong,mnlat,nlev)           :: global
      REAL, DIMENSION(:,:,:),   ALLOCATABLE       :: ddepo
      REAL, DIMENSION(:,:,:),   ALLOCATABLE       :: sddat
! temporary extra variable for checks
      REAL                                        :: missings
      REAL                                        :: missingbe7
      REAL                                        :: missingbe10
      REAL                                        :: missinglnox
      REAL                                        :: missingacnox
      REAL                                        :: aotau

      CHARACTER(LEN=40), DIMENSION(numflux), SAVE :: fnames
      CHARACTER(LEN=10), PARAMETER                :: version='601.1'
      CHARACTER(LEN=14)                           :: filename
      CHARACTER(LEN=80)                           :: cmessage

! missing cell mask for lightning and volcano emissions
      LOGICAL, DIMENSION(:,:,:), ALLOCATABLE       :: mask

      COMMON /sum3_d/global                  ! Used in SUM3D

!-----------------------------------------------------------------------

      ALLOCATE(cellswap(2,nclprc*nproc))
      ALLOCATE(ipos(3,nclprc))
      ALLOCATE(nbl(nlnpe,nlpe))
      ALLOCATE(nnn(nlnpe,nlpe,nlev))
      ALLOCATE(nph(nlnpe,nlpe,nlev))

      ALLOCATE(posold(4,nclprc))
      ALLOCATE(dja(ndj,nlnpe,nlpe,nlev))
      ALLOCATE(dj(ndj,nclprc))
      ALLOCATE(cloud(nlnpe,nlpe,nlev))
      ALLOCATE(seaicefr(nlnpe,nlpe))
      ALLOCATE(snowfr(nlnpe,nlpe))
      ALLOCATE(psurf(nlnpe,nlpe))
      ALLOCATE(tsurf(nlnpe,nlpe))
      ALLOCATE(tropdat(nlnpe,nlpe))
      ALLOCATE(landum(nlonpe,nlatpe))
      ALLOCATE(tl(nclprc))
      ALLOCATE(ql(nclprc))
      ALLOCATE(rhl(nclprc))
      ALLOCATE(m(nclprc))
      ALLOCATE(liq(nclprc))
      ALLOCATE(em(nc+2,nclprc)) ! includes L & Ac NOx
      ALLOCATE(dd(nc,nclprc))
      ALLOCATE(dw(nc,nclprc))
      ALLOCATE(mass(numchem,nlnpe,nlpe,nlev))
      ALLOCATE(m1(numchem,nlnpe,nlpe,nlev))
      ALLOCATE(ddepo(nc,nlnpe,nlpe))
      ALLOCATE(sddat(nc+1,nlnpe,nlpe))
      ALLOCATE(rc(nlonpe,nlatpe,ntype,nc))
      ALLOCATE(o3_stom_frac(nlonpe,nlatpe))
      ALLOCATE(rh(nlonpe,nlatpe,nmetlev))
      ALLOCATE(qs(nlonpe))

      ALLOCATE(terpem(nlnpe,nlpe))

      ALLOCATE(mask(nlnpe,nlpe,nlev))

! DEPENDS ON: cpustats
      CALL CPUSTATS('MAIN        ',1)
      WRITE(6,*) 'STOCHEM called at: ',day,month,year
! The UM may return a day equal to days of month + 1, so correct.
      IF (day >= (daym(month)+1.0-1.0e-6)) THEN
        day = day - daym(month)
        month = month + 1
        IF (month > 12) THEN
          month = 1
          year = year + 1
! DEPENDS ON: set_daym
          CALL SET_DAYM(daym,year)
          ysec = REAL(SUM(daym)*24*3600)
        END IF
        WRITE(6,*) 'Reset DAY on entry to:',day,month,year
      END IF

! Only allow start from dump on step 1
      IF (UMStepNo /= 1) StochstartDump=.false.

! Base from which to start counting cell numbers, PE 0 has 1->NCLPRC,
! PE 1 has NCLPRC+1->2*NCLPRC etc.
      cellbase=mype*nclprc

! Initialise model
      IF (mype == 0) WRITE(6,*) 'NSTEP=',nstep
      IF (nstep == 0) then
! DEPENDS ON: cpustats
        CALL CPUSTATS('INIT_STOCHEM',1)
        CALL INITIALISE_STOCHEM
        CALL KK_GCOM_INI(nclprc*(nc+2))
! DEPENDS ON: cpustats
        CALL CPUSTATS('INIT_STOCHEM',2)
      ELSE
! DEPENDS ON: secs
        time = SECS(day,month,year)
! DEPENDS ON: set_daym
        CALL SET_DAYM(daym,year)       ! in case its a new year
        ysec = REAL(SUM(daym)*daysec)  ! for leap years etc
        max_cells = 0
      END IF

! Set filetypes for pp output
      IF (lhourly) THEN
        filetype2='a'
      ELSE
        filetype2='m'
      END IF

      endtime=time+stochem_call_step
      WRITE(6,*) 'TIME,ENDTIME=',time,endtime

!
! Calculate relative humidity. Do outside of main loop below as
! temperature, pressure will not change
      DO l = 1, nmetlev
        DO k = 1, rowspe
! DEPENDS ON: qsat
          CALL QSAT(qs,t(:,k,l),p(:,k,l),nlonpe)
! Calculate relative humidity (fraction)
          DO i = 1, nlonpe
            rh(i,k,l) = q(i,k,l) / qs(i)
            rh(i,k,l) = MIN(rh(i,k,l),1.0)
          END DO
        END DO
      END DO

! Main advection timestep loop start
! DEPENDS ON: cpustats
      CALL CPUSTATS('STEP        ',1)

      DO
        IF (time >= endtime) EXIT

! cellflux = 0.0 now done in ADVFLUX
!       cellflux=0.

        WRITE(6,*) 'mype=',mype,' DAY: ',day,' MONTH: ',month,          &
     &    ' YEAR: ',year

! Dry Depositions
! DEPENDS ON: cpustats
        CALL CPUSTATS('DRYDEP      ',1)
! DEPENDS ON: surf_ddep_res
        CALL SURF_DDEP_RES(t0,p(:,:,0),rh(:,:,1),soilmc,so4_vd,gsf,     &
     &    stcon,t0tile,lai_ft,canwc,rc,o3_stom_frac)
!       CALL CPUSTATS('DRYDEP      ',2)

!       CALL CPUSTATS('DRYDEP      ',1)
! DEPENDS ON: drydep
        CALL DRYDEP(bl,gsf,ra,rq,rc,o3_stom_frac,ddepo,sddat)
! DEPENDS ON: cpustats
        CALL CPUSTATS('DRYDEP      ',2)

! DEPENDS ON: cpustats
        CALL CPUSTATS('SURFCALC    ',1)
! DEPENDS ON: surfcalc
        CALL SURFCALC(sice,sa,seaicefr,snowfr)
! DEPENDS ON: cpustats
        CALL CPUSTATS('SURFCALC    ',2)

        posold=pos

! Advection: reposition cells with/without time interpolation
! DEPENDS ON: cpustats
        CALL CPUSTATS('RKO4_SINGLE ',1)
! DEPENDS ON: rko4_single_time
        CALL RKO4_SINGLE_TIME(pos,u,v,w,bl,orog,seed2(mype,:),nfill,    &
     &    cellswap,cellbase,z_top_of_model,first_constant_r_rho_level)
! DEPENDS ON: cpustats
        CALL CPUSTATS('RKO4_SINGLE ',2)

! DEPENDS ON: cpustats
        CALL CPUSTATS('ADVFLUX     ',1)
! DEPENDS ON: advflux
        CALL ADVFLUX(pos,tropz,xx,flist,totflu,cellflux,nflux,posold,   &
     &    nfill,lnp)
! DEPENDS ON: cpustats
        CALL CPUSTATS('ADVFLUX     ',2)

! DEPENDS ON: cpustats
        CALL CPUSTATS('SWAP        ',1)
#ifdef OLDSWAP1
        CALL SWAP ! uses CELLSWAP to swap POS and XX
#else
        CALL EXCH_SWAP
#endif
! DEPENDS ON: cpustats
        CALL CPUSTATS('SWAP        ',2)

! DEPENDS ON: blrand
        CALL BLRAND(pos,bl,orog,z_top_of_model,                         &
     &    first_constant_r_rho_level,nfill,seed2(mype,:))

! DEPENDS ON: cpustats
        CALL CPUSTATS('CLMIX2      ',1)
! DEPENDS ON: clmix2
        CALL CLMIX2(pos,p,lnp,up,upd,orog,Stochem_advection_step,       &
     &    seed3(mype,:),nfill,z_top_of_model,first_constant_r_rho_level)
! DEPENDS ON: cpustats
        CALL CPUSTATS('CLMIX2      ',2)

! DEPENDS ON: cpustats
        CALL CPUSTATS('GETNNN      ',1)
! DEPENDS ON: getnnn
        CALL GETNNN(ipos,pos,nnn,nbl,bl,orog,nfill)
! DEPENDS ON: cpustats
        CALL CPUSTATS('GETNNN      ',2)

! DEPENDS ON: cpustats
        CALL CPUSTATS('MISS        ',1)
! DEPENDS ON: miss
        CALL MISS(missings,missingbe7,missingbe10,missinglnox,          &
     &    missingacnox,so2em,be7em,be10em,lnoxem,acnoxem,               &
     &    totso2em,totbe7em,totbe10em,totlnoxem,                        &
     &    totacnoxem,nnn)
! DEPENDS ON: cpustats
        CALL CPUSTATS('MISS        ',2)

! Establish photolysis rates at start, 1/3, 2/3 and end of Advection STE
! 10 day timescale
        aotau=Stochem_advection_step/strat_tau

! Change the ozone profile slowly, exp(aotau)=1- aotua for small values
        WHERE(nnn/=0) o3conc=o3conc*(1.0-aotau)+conc(i_o3,:,:,:)*aotau
        WHERE(nnn/=0) ch4conc=ch4conc*(1.0-aotau)+                      &
     &    conc(i_ch4,:,:,:)*aotau

! Use O3 climatology above tropopause for photolysis rates
! DEPENDS ON: cpustats
        CALL CPUSTATS('STRATMASK   ',1)
! DEPENDS ON: stratmask
        CALL STRATMASK(tropz,o3um,orog,o3conc,z_top_of_model,           &
     &    first_constant_r_rho_level)
! DEPENDS ON: cpustats
        CALL CPUSTATS('STRATMASK   ',2)

! DEPENDS ON: cpustats
        CALL CPUSTATS('CLCALC      ',1)
! DEPENDS ON: clcalc
        CALL CLCALC(cloud,dc,cc)
! DEPENDS ON: cpustats
        CALL CPUSTATS('CLCALC      ',2)

! DEPENDS ON: cpustats
        CALL CPUSTATS('JVALS       ',1)
! DEPENDS ON: jvals
        CALL JVALS(dj,dja,ipos,cloud,o3col,o3conc,land,snowfr,seaicefr, &
     &    time,t0,t,p,nfill)
! DEPENDS ON: cpustats
        CALL CPUSTATS('JVALS       ',2)

! Get temperature, humidity and cloud water for this cell by
! 3D interpolation
! DEPENDS ON: cpustats
        CALL CPUSTATS('TEMP        ',1)
! DEPENDS ON: temp
        CALL TEMP(pos,t,t0,tl,nfill)
! DEPENDS ON: cpustats
        CALL CPUSTATS('TEMP        ',2)
! DEPENDS ON: cpustats
        CALL CPUSTATS('WATER       ',1)
! DEPENDS ON: water
        CALL WATER(pos,q,clw,rh,ql,liq,rhl,nfill)
! DEPENDS ON: cpustats
        CALL CPUSTATS('WATER       ',2)

! Calculate molecular density of air in molecule/cm3
! DEPENDS ON: cpustats
        CALL CPUSTATS('MOLDENSE    ',1)
! DEPENDS ON: moldense
        CALL MOLDENSE(m,pos,lnp,tl,nfill)
! DEPENDS ON: cpustats
        CALL CPUSTATS('MOLDENSE    ',2)

! Get b.l. depth, Va, precipitation and cloud top by nearest neighbour
! DEPENDS ON: cpustats
        CALL CPUSTATS('EGET        ',1)
! DEPENDS ON: eget
        CALL EGET(estore,em,emiss,ipos,pos,bl,orog,nnn,nbl,month,       &
     &    time,so2em,missings,be7em,missingbe7,be10em,                  &
     &    missingbe10,lnoxem,missinglnox,acnoxem,missingacnox,nfill)
! DEPENDS ON: cpustats
        CALL CPUSTATS('EGET        ',2)

! DEPENDS ON: cpustats
        CALL CPUSTATS('DDCALC      ',1)
! DEPENDS ON: ddcalc
        CALL DDCALC(dstore,dd,ddepo,ipos,pos,bl,orog,nbl,nfill)
! DEPENDS ON: cpustats
        CALL CPUSTATS('DDCALC      ',2)

! DEPENDS ON: cpustats
        CALL CPUSTATS('DWCALC      ',1)
! DEPENDS ON: dwcalc
        CALL DWCALC(dw,adp,acp,cc,pos,nfill,lnp)
! DEPENDS ON: cpustats
        CALL CPUSTATS('DWCALC      ',2)

! Calculate stratospheric input every 3 hours
! DEPENDS ON: cpustats
        CALL CPUSTATS('STRATCALC2  ',1)
!       CALL STRATCALC2(xx,em,pos,ipos,nfill,trop,
!    &    o3_ls(:,:,1:nmetlev),stochem_advection_step,lnp)
! DEPENDS ON: stratcalc3
        CALL STRATCALC3(xx,m,em,pos,ipos,nfill,tropz,orog,o3um,         &
     &    stochem_advection_step)

! DEPENDS ON: cpustats
        CALL CPUSTATS('STRATCALC2  ',2)

! Swap cells between PEs for more efficient CPU usage
! DEPENDS ON: cpustats
        CALL CPUSTATS('POOL        ',1)
#ifdef OLDSWAP2
        CALL POOL  ! sets up CELLSWAP (passed in COMMON block)
#else
        CALL MAKE_SWAP_MAP
#endif
! DEPENDS ON: cpustats
        CALL CPUSTATS('POOL        ',2)
! DEPENDS ON: cpustats
        CALL CPUSTATS('SWAP2       ',1)
#ifdef OLDSWAP2
        CALL SWAP2 ! uses CELLSWAP to swap POS and XX etc
#else
        CALL BALANCE_SWAP(.TRUE.)
#endif
! DEPENDS ON: cpustats
        CALL CPUSTATS('SWAP2       ',2)

! Calculate new concentrations
! DEPENDS ON: cpustats
        CALL CPUSTATS('DERIV       ',1)
! DEPENDS ON: deriv
        CALL DERIV(nflux,nfill,flist,stochem_advection_step,time,tl,    &
     &    ql,rhl,m,liq,dj,dd,dw,em,xx,totflu,cellflux)

! DEPENDS ON: cpustats
        CALL CPUSTATS('DERIV       ',2)

! Swap all the cells back to their original processors
#ifdef OLDSWAP2
! exchanges elements 1 and 2 of CELLSWAP
        cellswap=cshift(cellswap,shift=1,dim=1)
! DEPENDS ON: cpustats
        CALL CPUSTATS('SWAP        ',1)
        CALL SWAP ! uses CELLSWAP to swap POS and XX
! DEPENDS ON: cpustats
        CALL CPUSTATS('SWAP        ',2)
#else
! DEPENDS ON: cpustats
        CALL CPUSTATS('SWAP2       ',1)
        CALL BALANCE_SWAP(.FALSE.)
! DEPENDS ON: cpustats
        CALL CPUSTATS('SWAP2       ',2)
#endif

! DEPENDS ON: follow
        IF(lfol) CALL FOLLOW(pos,xx,cellflux,clist,flist,follist,       &
     &    cellno,nchem,nflux,nfill,time,day,month,year)

! Generate ipos array
! DEPENDS ON: cpustats
        CALL CPUSTATS('AINDEX      ',1)
! DEPENDS ON: aindex
        CALL AINDEX(nfill,pos,ipos)
! DEPENDS ON: cpustats
        CALL CPUSTATS('AINDEX      ',2)

! Add up 3D fluxes. First, identify which 3D fluxes are required.
        n3dflux = 0
        DO n = 1, nflux
          IF (flist(2,n) > 0 .AND. flist(1,n) < 1701) THEN
            n3dflux = n3dflux + 1
            idxflux(n3dflux) = flist(2,n)
          END IF
        END DO

        DO j=1,nfill
          i=ipos(1,j)
          k=ipos(2,j)
          l=ipos(3,j)
!CDIR NODEP
          DO n=1,n3dflux
            tflux(idxflux(n),i-lndat+1,k-ltdat+1,l) =                   &
     &        tflux(idxflux(n),i-lndat+1,k-ltdat+1,l) +                 &
     &        cellflux(idxflux(n),j)
          END DO
        END DO
!
! AOT 40. Convert units from ppbv-s to ppbv-h.
        DO n=1,nflux
          k2 = flist(2,n)
          IF (flist(1,n) == 1799 .AND. k2 > 0) THEN
            DO j=1,nfill
              i=ipos(1,j)-lndat+1
              k=ipos(2,j)-ltdat+1
              l=ipos(3,j)
              IF (nnn(i,k,l) > 0) THEN
                tflux(k2,i,k,l) = tflux(k2,i,k,l) +                     &
     &            (cellflux(k2,j) / REAL(nnn(i,k,l))) *                 &
     &            stochem_advection_step / 3600.0
              END IF
            END DO
          END IF
        END DO
!
! Surface deposition velocities (1/rc), o3 stomatal fraction
!CDIR NODEP
        DO n=1,nflux
          IF (flist(1,n) > 1700 .AND. flist(1,n) <= (1700+nc+1)) THEN
            IF (flist(2,n) > 0) THEN
              DO k = 1, nlpe
                DO i = 1, nlnpe
                  tflux(flist(2,n),i,k,1) = tflux(flist(2,n),i,k,1) +   &
     &              sddat(flist(1,n)-1700,i,k)
                  totflu(n) = totflu(n) + sddat(flist(1,n)-1700,i,k)
                END DO
              END DO
            END IF
! Add on lightning flash counts
          ELSE IF (flist(1,n) > 1795 .AND. flist(1,n) < 1798) THEN
            IF (flist(2,n) > 0) THEN
              DO k = 1, nlpe
                DO i = 1, nlnpe
                  tflux(flist(2,n),i,k,1) = tflux(flist(2,n),i,k,1) +   &
     &              nflashes(i,k,flist(1,n)-1795)
                  totflu(n) = totflu(n) + nflashes(i,k,flist(1,n)-1795)
                END DO
              END DO
            END IF
          END IF
        END DO

! Place emissions/depositions into storage if no cells are present
! or reset storage
! DEPENDS ON: cpustats
        CALL CPUSTATS('STORE       ',1)
! DEPENDS ON: store
        CALL STORE(estore,emiss,nnn,nbl,time,month)
! DEPENDS ON: cpustats
        CALL CPUSTATS('STORE       ',2)
! DEPENDS ON: cpustats
        CALL CPUSTATS('STORED      ',1)
! DEPENDS ON: stored
        CALL STORED(dstore,ddepo,nbl)
! DEPENDS ON: cpustats
        CALL CPUSTATS('STORED      ',2)

! Map Lagrangian concs. onto Eulerian grid and do mixing
! DEPENDS ON: cpustats
        CALL CPUSTATS('PLOT        ',1)
! DEPENDS ON: plot
        CALL PLOT(conc,xx,ipos,nnn,mix_fact1,mix_fact2,nfill)
! DEPENDS ON: cpustats
        CALL CPUSTATS('PLOT        ',2)

! Collect statistics
! DEPENDS ON: cpustats
        CALL CPUSTATS('STATS       ',1)
! DEPENDS ON: stats
        CALL STATS(mconc,sdconc,nm,np,nph,nstat,nnn,conc,clist,nchem)
! DEPENDS ON: cpustats
        CALL CPUSTATS('STATS       ',2)
! DEPENDS ON: cpustats
        CALL CPUSTATS('MCALC       ',1)
! DEPENDS ON: mcalc
        CALL MCALC(mass,totmas,xx,clist,ipos,pos,tropz,nchem,nfill)
! DEPENDS ON: cpustats
        CALL CPUSTATS('MCALC       ',2)
        navg = navg + 1
        totavg = totavg + totmas

        time = time + stochem_advection_step

! Calculate day and month. If month changes, get new emissions.
        day = day + stochem_advection_step / daysec
        period = period + stochem_advection_step / daysec
        IF (day >= (daym(month)+1.0-1.0e-6)) THEN ! n.b. tolerance

! Produce monthly output
          IF (.NOT.lhourly) THEN
            WRITE(6,*) 'Monthly_Output Called: ',day,month,year
! DEPENDS ON: cpustats
            CALL CPUSTATS('MONTHLY_OUTP',1)
            CALL MONTHLY_OUTPUT
! DEPENDS ON: cpustats
            CALL CPUSTATS('MONTHLY_OUTP',2)
! DEPENDS ON: cpustats
            CALL CPUSTATS('ZERO_STATS  ',1)
            CALL ZERO_STATS
! DEPENDS ON: cpustats
            CALL CPUSTATS('ZERO_STATS  ',2)
          END IF
          day = day - daym(month)
          daym0 = day
          month = month + 1
          IF (month > 12) THEN
            month = 1
            year = year + 1
! DEPENDS ON: set_daym
            CALL SET_DAYM(daym,year)
            ysec = REAL(SUM(daym)*daysec)
          END IF

! Initialise new monthly station file
          IF (mype == 0) THEN
            IF (lstat) THEN
! DEPENDS ON: findname
              CALL FINDNAME('s','m','c',0,0,filename,month,year)
              OPEN(59,FILE=filename,STATUS='REPLACE')
              WRITE(6,*) 'New Station file created: ',filename
              WRITE(59,*) nstation
              DO i=1,nstation
                WRITE(59,*) station_names(i)
              END DO
              WRITE(59,*) nchem
              DO i=1,nchem
                WRITE(59,*) cnames(clist(i))
              END DO
              CLOSE(59)
            END IF

! Open file for cell output
            IF (lfol) THEN
! DEPENDS ON: findname
              CALL FINDNAME('f','m','c',0,0,filename,month,year)
              OPEN(61,FILE=filename,STATUS='REPLACE')
              WRITE(6,*) 'New cell following file created: ',filename
              WRITE(61,*) numfol
              CLOSE(61)
            END IF
          END IF

          IF (mype==0) THEN
            WRITE(6,*)'New Month, call EMCALC & lightread & o3read'
            WRITE(6,'(a5,2(i2,a1),i4)')                                 &
     &        'Day: ',INT(day),'-',month,'-',year
          END IF

! Call appropriate emissions routine.
          SELECT CASE(scenario)
          CASE('fi')
! DEPENDS ON: emcalc_a1fi
            CALL EMCALC_A1FI(emiss,month,year)
          CASE('ab')
! DEPENDS ON: emcalc_a1b
            CALL EMCALC_A1B(emiss,month,year)
          CASE('a2')
! DEPENDS ON: emcalc_a2
            CALL EMCALC_A2(emiss,month,year)
          CASE('b1')
! DEPENDS ON: emcalc_b1
            CALL EMCALC_B1(emiss,month,year)
          CASE('b2')
! DEPENDS ON: emcalc_b2
            CALL EMCALC_B2(emiss,month,year)
          CASE('bu')
! DEPENDS ON: emcalc_cle
            CALL EMCALC_CLE(emiss,month,year)
          CASE('mf')
! DEPENDS ON: emcalc_mfr
            CALL EMCALC_MFR(emiss,month,year)
          CASE('A2')
! DEPENDS ON: emcalc_a2_
            CALL EMCALC_A2_(emiss,month,year)
          CASE('pi')
! DEPENDS ON: emcalc_pi
            CALL EMCALC_PI(emiss,month,year)
          CASE DEFAULT
            cmessage='NO EMISSION SCENARIO SELECTED'
! DEPENDS ON: ereport
            CALL EREPORT('STOCHEM',1,cmessage)
          END SELECT

! Read in 72x36x5 Li & Shine O3 climatology for the month
! DEPENDS ON: cpustats
          CALL CPUSTATS('READO3      ',1)
! DEPENDS ON: reado3
          CALL READO3(o3_ls(:,:,1:nmetlev),month,z_top_of_model,        &
     &      first_constant_r_rho_level)
! DEPENDS ON: cpustats
          CALL CPUSTATS('READO3      ',2)

! Calculate a column O3 value (0-model top) : o3col (Dobson units)
! DEPENDS ON: cpustats
          CALL CPUSTATS('O3COLUMN    ',1)
! DEPENDS ON: o3column
          CALL O3COLUMN(o3_ls(:,:,1:nmetlev),o3col)
! DEPENDS ON: cpustats
          CALL CPUSTATS('O3COLUMN    ',2)
! DEPENDS ON: cpustats
          CALL CPUSTATS('LIGHTREAD   ',1)
! DEPENDS ON: lightread
          CALL LIGHTREAD(acnoxem,so2em,be7em,be10em,month,year,         &
     &      totacnoxem,totso2em,totbe7em,totbe10em,p,lnp)
! DEPENDS ON: cpustats
          CALL CPUSTATS('LIGHTREAD   ',2)
        END IF  ! Month change

! Calculate Lightning NOx emissions
        IF (lightningon) THEN
! DEPENDS ON: lightnox
          CALL LIGHTNOX(cc_base_level,cc_top_level,t,acr,orog,          &
     &      land_fraction,totlnoxem,lnoxem,nflashes,z_top_of_model,     &
     &      first_constant_r_rho_level)
        ELSE
          lnoxem = 0.0
          totlnoxem = 0.0
          nflashes = 0.0
        END IF

! Calculate isoprene and monoterpene emissions
! DEPENDS ON: emcalc_natural_voc
        CALL EMCALC_NATURAL_VOC(month,day,totpar,dirpar,gsf,t0tile,     &
     &    lai_ft,emiss(i_c5h8,:,:),terpem)

        nstep = nstep + 1
        IF (mype==0) WRITE(6,*) 'NSTEP: ',nstep

! Station data
        IF (lstat) THEN
! DEPENDS ON: station
          CALL STATION(xx,pos,cellno,bl,orog,stlon,stlat,nstation,day,  &
     &      month,year,clist,nchem,nfill)
        END IF

! Map Lagrangian concs. onto Eulerian grid and do mixing
! DEPENDS ON: cpustats
        CALL CPUSTATS('PLOT        ',1)
! DEPENDS ON: plot
        CALL PLOT(conc,xx,ipos,nnn,0.,0.,nfill)
! DEPENDS ON: cpustats
        CALL CPUSTATS('PLOT        ',2)

! Calculate the maximum number of cells carried on a processor.
        WRITE(6,*) 'Max cells on PE:',mype,max_cells(mype)
        IF(mype == 0) THEN
          WRITE(6,*) 'DAY: ',day,'Maximum number of cells: ',           &
     &      MAXVAL(max_cells),' PE: ',MAXLOC(max_cells)
          max_cells=0
        END IF

! Main advection timestep loop finish.
      END DO

! DEPENDS ON: cpustats
      CALL CPUSTATS('STEP        ',2)

! Produce hourly output
      IF (lhourly .AND. MOD(time,hoursec) < 1.0) THEN
        CALL MONTHLY_OUTPUT
        CALL ZERO_STATS
        daym0=day
      END IF
      IF (StochDumpRequest) THEN
! DEPENDS ON: cpustats
        CALL CPUSTATS('DUMP        ',1)
        CALL DUMP
! DEPENDS ON: cpustats
        CALL CPUSTATS('DUMP        ',2)
      END IF

! Calculate 3D fields of O3,CH4
! DEPENDS ON: data2met3d
      CALL DATA2MET3D(o3conc,o3mmr,mo3)
! DEPENDS ON: data2met3d
      CALL DATA2MET3D(ch4conc,ch4mmr,mch4)

! DEPENDS ON: cpustats
      CALL CPUSTATS('MAIN        ',2)
! DEPENDS ON: cpustats
      CALL CPUSTATS('MAIN        ',3)

      DEALLOCATE(cellswap)
      DEALLOCATE(ipos)
      DEALLOCATE(nbl)
      DEALLOCATE(nnn)
      DEALLOCATE(nph)

      DEALLOCATE(posold)
      DEALLOCATE(dja)
      DEALLOCATE(dj)
      DEALLOCATE(cloud)
      DEALLOCATE(seaicefr)
      DEALLOCATE(snowfr)
      DEALLOCATE(psurf)
      DEALLOCATE(tsurf)
      DEALLOCATE(tropdat)
      DEALLOCATE(landum)
      DEALLOCATE(tl)
!      DEALLOCATE(aph)
      DEALLOCATE(ql)
      DEALLOCATE(m)
      DEALLOCATE(rhl)
      DEALLOCATE(liq)
      DEALLOCATE(em)
      DEALLOCATE(dd)
      DEALLOCATE(dw)
      DEALLOCATE(mass)
      DEALLOCATE(m1)
      DEALLOCATE(ddepo)
      DEALLOCATE(sddat)
      DEALLOCATE(mask)
      DEALLOCATE(terpem)
      DEALLOCATE(rc)
      DEALLOCATE(o3_stom_frac)
      DEALLOCATE(rh)
      DEALLOCATE(qs)

!#######################################################################
      CONTAINS
!#######################################################################

!#######################################################################

      SUBROUTINE BALANCE_SWAP(forward)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Swap cells forth and back for load balancing
!-                         Internal procedure to STOCHEM
!-
!-   Inputs  :
!-   Outputs :
!-   Controls:
!-
! Current Code Owner: M.G. Sanderson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   06/04/05   Created. R. Johanni
!
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
      LOGICAL :: forward

! cellno and pos are not needed in deriv, so no need
! to swap them here

      CALL SWAP_ARRAY(xx,SIZE(xx,1),forward)
      CALL SWAP_ARRAY(cellflux,SIZE(cellflux,1),forward)

      IF (forward) THEN
        CALL SWAP_ARRAY(tl,1,forward)
        CALL SWAP_ARRAY(ql,1,forward)
        CALL SWAP_ARRAY(liq,1,forward)
        CALL SWAP_ARRAY(m,1,forward)
        CALL SWAP_ARRAY(em,SIZE(em,1),forward)
        CALL SWAP_ARRAY(dd,SIZE(dd,1),forward)
        CALL SWAP_ARRAY(dw,SIZE(dw,1),forward)
        CALL SWAP_ARRAY(dj,SIZE(dj,1),forward)
        nfill = nc_target(mype)
      ELSE
        nfill = nc_actual(mype)
      ENDIF

      END SUBROUTINE BALANCE_SWAP
      SUBROUTINE DUMP
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Dumps program variables to file.
!-                          Internal procedure to STOCHEM
!-
!-   Inputs  : Indirectly- XX,POS,ESTORE,DSTORE,TFLUX,TOTFLU,FLIST,
!-             FNAMES,CONC,MCONC,SDCONC,O3CONC,NM,NP,M0
!-             TIME,YEAR,MONTH,DAY,SEED2,SEED3,NFLUX
!-   Outputs : none
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.0    17/01/95  Created.  W.J. Collins
!  4.4    12/10/98  Added O3CONC to dump. D.S. Stevenson
!  4.5    25/03/99  Added UM Step No. to dump. C.E. Johnson
!  4.5    11/05/99  Now uses UM-style file names. C.E. Johnson
!  5.5    23/09/04  Changed to save seeds as integer arrays plus
!                   other variables added to list. M.G. Sanderson.
!  6.1    20/10/04  No change.
!

!VVV  V5.1  DUMP  2/8/01 - B deleted
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------

      INTEGER :: info
      INTEGER :: imsg=1
      INTEGER :: j
      INTEGER :: i
      CHARACTER(LEN=14) :: filename

! Open new file on PE 0 - write cell information
      info=GC_SHM_PUT
      CALL GC_SSYNC(nproc,info)
! DEPENDS ON: findname
      CALL FINDNAME('d','a','c',0,0,filename,                           &
     &              month,year,INT((day-INT(day))*24),INT(day))
      IF(mype==0) THEN
        WRITE(6,*) 'Writing ',filename,' at ',time,' day ',day,         &
     &    ' month ',month,' UM Step No. ',UMStepNo
        OPEN(51,FILE=filename,STATUS='UNKNOWN',                         &
     &    FORM='UNFORMATTED')
        WRITE(51) time,year,month,day
        WRITE(51) nproc,nlong,mnlat,nlev,ncell,nc
        WRITE(51) numflux,num3dflux,num3dflux_dim,numchem,outnames,     &
     &    fluxnames
        DO j=1,nfill
          WRITE(51) mype,j,pos(:,j),xx(:,j),cellno(j)
        END DO
        CLOSE(51)
      END IF

! Append to dump file from other PEs in turn - write cell info
      DO i=1,nproc-1
        info=GC_SHM_PUT
        CALL GC_SSYNC(nproc,info)
        IF(mype==i) THEN
          OPEN(51,FILE=filename,STATUS='OLD',                           &
     &      FORM='UNFORMATTED',POSITION='APPEND')
          DO j=1,nfill
            WRITE(51) mype,j,pos(:,j),xx(:,j),cellno(j)
          END DO
          CLOSE(51)
        END IF
      END DO

! Sum TOTAVG and TOTFLU over all processors (also SEEDs)
      info=GC_SHM_PUT
      CALL GC_SSYNC(nproc,info)
      CALL GC_RSUM(2*nc,nproc,info,totavg)
      CALL GC_RSUM(2*nc,nproc,info,totm0)
      CALL GC_RSUM(numflux,nproc,info,totflu)
! Since SEEDn(mype) is only non-zero element. Sum is equivalent to
! broadcst over PEs.
      CALL GC_ISUM(maxproc*ransize,nproc,info,seed2)
      CALL GC_ISUM(maxproc*ransize,nproc,info,seed3)
      info=GC_SHM_PUT
      CALL GC_SSYNC(nproc,info)
      IF(mype/=0) THEN
! Need to keep sum over all PEs unchanged, so set all others to zero
        totm0 = 0.0
        totavg = 0.0
        totflu = 0.0
      END IF

! Append Eulerian data to dump file, sent from PE0
      IF(mype==0) OPEN(51,FILE=FILENAME,STATUS='OLD',                   &
     &    FORM='UNFORMATTED',POSITION='APPEND')
      DO i=1,nc
! DEPENDS ON: sum2d
        CALL SUM2D(estore(i,:,:))
        IF(mype==0) WRITE(51) global(:,:,1)
! DEPENDS ON: sum2d
        CALL SUM2D(dstore(i,:,:))
        IF(mype==0) WRITE(51) global(:,:,1)
! DEPENDS ON: sum3d
        CALL SUM3D(conc(i,:,:,:))
        IF(mype==0) WRITE(51) global
      END DO
      DO i=1,numchem
! DEPENDS ON: sum3d
        CALL SUM3D(mconc(i,:,:,:))
        IF(mype==0) WRITE(51) global
! DEPENDS ON: sum3d
        CALL SUM3D(sdconc(i,:,:,:))
        IF(mype==0) WRITE(51) global
! DEPENDS ON: sum3d
        CALL SUM3D(m0(i,:,:,:))
        IF(mype==0) WRITE(51) global
      END DO
      DO i=1,num3dflux_dim
! DEPENDS ON: sum3d
        CALL SUM3D(tflux(i,:,:,:))
        IF(mype==0) WRITE(51) global
      END DO
! DEPENDS ON: sum3d
      CALL SUM3D(o3conc)
      IF(mype==0) WRITE(51) global
! DEPENDS ON: sum3d
      CALL SUM3D(REAL(nm))
      IF(mype==0) WRITE(51) global
! DEPENDS ON: sum3d
      CALL SUM3D(REAL(np))
      IF(mype==0) WRITE(51) global
      IF(mype==0) CLOSE(51)

! Append all other data to dump file, sent from PE 0
      IF(mype==0) THEN
        OPEN(51,FILE=filename,STATUS='OLD',                             &
     &    FORM='UNFORMATTED',POSITION='APPEND')
        WRITE(51) totflu
        WRITE(51) flist,fnames
        WRITE(51) stash_code
        WRITE(51) seed2
        WRITE(51) seed3
        WRITE(51) nflux,nstat
        WRITE(51) totm0
        WRITE(51) totavg
        WRITE(51) navg
        WRITE(51) daym0
        WRITE(51) Annual_Soil_NOx
        WRITE(51) Annual_Land_NVOC
        WRITE(51) Annual_Ocean_NVOC
        WRITE(51) Annual_Isop
        WRITE(51) Annual_DMS
        CLOSE(51)
      END IF

! Set the values of SEEDs back to zero that aren't for this PE.
      DO i = 0, maxproc-1
        IF (i /= mype) THEN
          seed2(i,:) = 0
          seed3(i,:) = 0
        END IF
      END DO

      END SUBROUTINE DUMP
!######################################################################

!#######################################################################

      SUBROUTINE EXCH_SWAP
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Swaps cell information between processors.
!-                         Internal procedure to STOCHEM.
!-
!-   Inputs  : also CELLSWAP, POS and XX
!-   Outputs : also POS and XX
!-   Controls:
!-
!
! Current Code Owner: M.G. Sanderson
!-----------------------------------------------------------------------
      USE MPL, ONLY :     &
          MPL_INTEGER,    &
          MPL_REAL,       &
          MPL_STATUS_SIZE

      IMPLICIT NONE

      INTEGER :: nswap, swapinfo(2,nclprc), num, npe, n, ier
      INTEGER :: nreqs, ireq(4*nproc), istat(MPL_STATUS_SIZE,4*nproc)
      INTEGER :: nsend(0:nproc-1), nrecv(0:nproc-1), n_items_recv
      INTEGER :: nx, np, nc
      INTEGER :: MY_COMM

      REAL, ALLOCATABLE :: tmp_x(:,:), tmp_p(:,:), tmp_c(:,:)
      INTEGER, ALLOCATABLE :: tmp_n(:)

! Get communicator for use in MPL comms
      CALL gc_get_communicator(MY_COMM,ier)

      nx = SIZE(xx,1)
      np = SIZE(pos,1)
      nc = SIZE(cellflux,1)

! Sort out which cells have to be swapped
! Have to remember to subtract 1 from CELLSWAP each time to get
! correct PE.

      nswap = 0
      DO n=1,nfill
        IF (cellswap(1,cellbase+n) > 0) THEN
          nswap = nswap+1
          swapinfo(1,nswap) = n
          swapinfo(2,nswap) = cellswap(2,cellbase+n)-1
        END IF
      END DO

! Count how many cells go to every PE

      nsend(:) = 0
      DO n=1,nswap
        nsend(swapinfo(2,n)) = nsend(swapinfo(2,n))+1
      END DO

      WRITE(6,*) 'EXCH_SWAP nsend: ',nsend

! Do an all-to-all comunication in order to see how many cells we
! receive

      CALL MPL_ALLTOALL(nsend,1,MPL_INTEGER,nrecv,1,MPL_INTEGER,        &
     &                  MY_COMM,ier)

      WRITE(6,*) 'EXCH_SWAP nrecv: ',nrecv

      n_items_recv = SUM(nrecv)

! Check if there is enough space allocated
      IF (nfill+n_items_recv > nclprc) THEN
        cmessage='STOCHEM Routine SWAP allocated too many cells '
        WRITE(6,*) cmessage
        WRITE(6,*) 'NFILL: ',NFILL,' N_ITEMS_RECV: ',N_ITEMS_RECV
        WRITE(6,*) 'NCLPRC: ',NCLPRC
! DEPENDS ON: ereport
        CALL EREPORT('SWAP',1,cmessage)
        CALL MPL_Abort(MY_COMM,1,ier)
      END IF

! Allocate temporary arrays for cells to be swapped out

      IF (nswap > 0) THEN
        ALLOCATE(tmp_x(nx,nswap))
        ALLOCATE(tmp_p(np,nswap))
        ALLOCATE(tmp_c(nc,nswap))
        ALLOCATE(tmp_n(nswap))
      END IF

! Write the cells to be swapped into the temporary arrays
! ordered by the PE where they go to

      num = 0
      DO npe=0,nproc-1
        IF (nsend(npe) > 0) THEN
          DO n=1,nswap
            IF (swapinfo(2,n) == npe) THEN
              num = num+1
              tmp_x(:,num) = xx(:,swapinfo(1,n))
              tmp_p(:,num) = pos(:,swapinfo(1,n))
              tmp_c(:,num) = cellflux(:,swapinfo(1,n))
              tmp_n(num)   = cellno(swapinfo(1,n))
            END IF
          END DO
        END IF
      END DO

! Safety check
      IF (num /= nswap) THEN
        WRITE(6,*) 'SWAP: internal error 1 num,nswap: ',num,nswap
        CALL MPL_ABORT(MY_COMM,1,ier)
      END IF

! Send the cells

      num = 0
      nreqs = 0 ! number of requests
      DO npe=0,nproc-1
        IF (nsend(npe) > 0) THEN
          CALL MPL_Isend(tmp_x(1,num+1),nx*nsend(npe),MPL_REAL,         &
     &                   npe,1,MY_COMM,ireq(nreqs+1),ier)
          CALL MPL_Isend(tmp_p(1,num+1),np*nsend(npe),MPL_REAL,         &
     &                   npe,2,MY_COMM,ireq(nreqs+2),ier)
          CALL MPL_Isend(tmp_c(1,num+1),nc*nsend(npe),MPL_REAL,         &
     &                   npe,3,MY_COMM,ireq(nreqs+3),ier)
          CALL MPL_Isend(tmp_n(num+1),nsend(npe),MPL_INTEGER,           &
     &                   npe,4,MY_COMM,ireq(nreqs+4),ier)
          nreqs = nreqs+4
          num = num+nsend(npe)
        END IF
      END DO

! Safety check
      IF (num /= nswap) THEN
        WRITE(6,*) 'SWAP: internal error 2 num,nswap: ',num,nswap
        CALL MPL_Abort(MY_COMM,1,ier)
      END IF

! Receive the cells

      DO npe=0,nproc-1
        IF (nrecv(npe) > 0) THEN
          CALL MPL_Recv(xx      (1,nfill+1),nx*nrecv(npe),MPL_REAL,     &
     &                  npe,1,MY_COMM,istat,ier)
          CALL MPL_Recv(pos     (1,nfill+1),np*nrecv(npe),MPL_REAL,     &
     &                  npe,2,MY_COMM,istat,ier)
          CALL MPL_Recv(cellflux(1,nfill+1),nc*nrecv(npe),MPL_REAL,     &
     &                  npe,3,MY_COMM,istat,ier)
          CALL MPL_Recv(cellno(nfill+1),nrecv(npe),MPL_INTEGER,         &
     &                  npe,4,MY_COMM,istat,ier)
          nfill = nfill+nrecv(npe)
        END IF
      END DO

! Wait for all communictation to complete

      IF (nreqs>0) CALL MPL_Waitall(nreqs,ireq,istat,ier)

! Deallocate temporary arrays

      IF (nswap>0) DEALLOCATE(tmp_x,tmp_p,tmp_c,tmp_n)

! Fill in holes using cells at end.
      DO n=nswap,1,-1
        num=swapinfo(1,n)
        xx(:,num)=xx(:,nfill)
        pos(:,num)=pos(:,nfill)
        cellflux(:,num)=cellflux(:,nfill)
        cellno(num)=cellno(nfill)
        nfill=nfill-1
      END DO

      END SUBROUTINE EXCH_SWAP
      SUBROUTINE GETLIS
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Get lists of fluxes and chemical species to
!-                         output. Internal procedure to STOCHEM.
!-
!-   Inputs  :
!-   Outputs : FLIST,FNAMES,NFLUX,CLIST,NCHEM,STLON,STLAT
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.0    30/06/94  Created.  W.J. Collins
!  5.1    27/10/00  Opens station and parcel following files as
!                   needed. C.E. Johnson
!  5.5    09/03/04  Now calculates string positions relative to
!                   len_flx_str parameter to avoid use of magic nos.
!                   M.G. Sanderson
!  6.2    01/03/06  Reads stash code from output strings.
!                                         M.G. Sanderson
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!----------------------------------------------------------------------

      INTEGER                               :: i,j,ierr,n3dflux
      INTEGER                               :: imsg=1
      INTEGER                               :: p1 = 41
      INTEGER                               :: p2 = 44
      INTEGER                               :: p3 = 47
      INTEGER                               :: p4 = 51
      INTEGER                               :: p5 = len_flx_str - 1
      CHARACTER(LEN=60)                     :: cstring
      CHARACTER*80                          :: cmessage
      CHARACTER(LEN=20), DIMENSION(numstat) :: stnames
      CHARACTER(LEN=14)                     :: filename
      LOGICAL                               :: lfile

      flist=0
      clist=0
      fnames=' '
      nchem=0
      ierr=0

! Construct CLIST using the OUTNAMES array
      DO j=1,numchem
        cstring = outnames(j)
        i=1
        DO
          IF (i > nc .OR. cstring(1:8) == cnames(i)) EXIT
          i=i+1
        END DO
        IF (i <= nc) THEN
          nchem = nchem + 1
          clist(nchem) = i
          IF (nchem >= numchem) EXIT
        END IF
      END DO

! Construct the fnames and flist arrays from fluxnames
! Assume reaction number and 3D output flag are stored at end
! of fluxnames string in format (I4,I2). p1, p2 and p3 are
! positions in string of these numbers.
      nflux = 0
      n3dflux = 0
      ierr = 0
      DO j=1,numflux
        cstring = fluxnames(j)
        IF(cstring(1:1)/='*') THEN
          nflux = nflux + 1
          READ(cstring,'(a40)') fnames(nflux)
          READ(cstring(p1:p2),'(i4)') flist(1,nflux)
          READ(cstring(p3:p4),'(i5)') stash_code(nflux)
          READ(cstring(p5:len_flx_str),'(i2)') flist(2,nflux)
          IF (flist(2,nflux) > 0) THEN
            n3dflux = n3dflux + 1
            IF (n3dflux > num3dflux) THEN
              WRITE(6,*) '**** ERROR GETLIS: N3DFLUX > NUM3DFLUX'
              WRITE(6,*) n3dflux,' > ',num3dflux
              cmessage='**** ERROR GETLIS: N3DFLUX > NUM3DFLUX'
! DEPENDS ON: ereport
              CALL EREPORT('GETLIS',1,cmessage)
            END IF
            flist(2,nflux)=n3dflux
          END IF
          IF (nflux > numflux) THEN
            cmessage='**** ERROR GETLIS: NFLUX > NUMFLUX'
            WRITE(6,*) cmessage
            WRITE(6,*) NFLUX,' > ',numflux
! DEPENDS ON: ereport
            CALL EREPORT('GETLIS',1,cmessage)
          END IF
        END IF
      END DO

! Construct the stnames, stlon, and stlat arrays from the
! Station_Names array.
      nstation=0
      ierr=0
      DO J=1,numstat
        cstring=station_names(j)
        IF (cstring(1:1)/='*') THEN
          nstation=nstation+1
          READ(cstring,'(a20)') stnames(nstation)
          READ(cstring(21:30),'(f10.1)') stlon(nstation)
          READ(cstring(31:40),'(f10.1)') stlat(nstation)
          IF (nstation >= numstat) EXIT
        END IF
      END DO

! Open new station file if needed
      IF (mype==0) THEN
        IF (LSTAT) THEN
! DEPENDS ON: findname
          CALL FINDNAME('s','m','c',0,0,FILENAME,MONTH,YEAR)
          INQUIRE(FILE=FILENAME,EXIST=LFILE)
          IF (.NOT. LFILE) THEN
            OPEN(59,FILE=FILENAME,STATUS='NEW')
            WRITE(6,*) 'GETLIS: New Station file created: ',FILENAME
            WRITE(59,*) NSTATION
            DO I=1,NSTATION
              WRITE(59,*) STNAMES(I)
            END DO
            WRITE(59,*) NCHEM
            DO I=1,NCHEM
              WRITE(59,*) CNAMES(CLIST(I))
            END DO
            CLOSE(59)
          END IF
        END IF
      END IF

! Construct the FOLLIST array
      FOLLIST=0
      DO I=1,NUMFOL
        FOLLIST(I)=cell_follow(i)
      END DO

! Open new cell following file if needed
      IF (mype==0) THEN
        IF (LFOL) THEN
! DEPENDS ON: findname
          CALL FINDNAME('f','m','c',0,0,FILENAME,MONTH,YEAR)
          INQUIRE(FILE=FILENAME,EXIST=LFILE)
          IF(.NOT. LFILE) THEN
            OPEN(61,FILE=FILENAME,STATUS='NEW')
            WRITE(6,*) 'GETLIS: New cell following file created: ',     &
     &        FILENAME
            WRITE(61,*) NUMFOL,' Cells followed'
            CLOSE(61)
          END IF
        END IF
      END IF

      END SUBROUTINE GETLIS
!#######################################################################
      SUBROUTINE INITIALISE_STOCHEM
!-----------------------------------------------------------------------
!
!     Initialise STOCHEM
!     Called this on first step of run.
!     Internal procedure to STOCHEM.
!
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.5    08/12/98  Created.  D.S. Stevenson
!  5.3    17/09/01  New Dynamics version. C.E. Johnson
!  5.5    13/01/04  Magic numbers removed. M.G. Sanderson
!  5.5    16/07/04  Now uses integer random number seeds, initialised
!                   using NEC routine FRANVI. M.G. Sanderson
!  5.5    02/08/04  Added scenario 'ab' to EMCALC list. M.G. Sanderson
!  6.2    01/03/06  Calculates initial random number seeds using
!                   routine FRANVI.  M.G. Sanderson
!
!
!VVV  V5.2 INITIALISE_STOCHEM  23/XI/00 DC,UP etc treated
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER :: info
      INTEGER :: imsg=1
      INTEGER :: i
      INTEGER :: j
      INTEGER :: mm
      INTEGER, DIMENSION(ransize) :: system_seed

      REAL :: ppb_vmr_o3 = 5.0e-9    ! ppbv for intervention
      REAL :: ppb_vmr_ch4 = 1.2e-6   ! ppbv for intervention
      REAL :: ppb_vmr = 1.0e-9       ! 1 ppbv in vmr
      REAL :: o3conc_init = 30.0e-09 ! Initial [O3] for photolysis in vm
      REAL :: ch4conc_init = 1.7e-6  ! Initial [CH4] in vmr
      REAL, DIMENSION(nlong,mnlat) :: temp_emiss

      CHARACTER(LEN=2) :: smonth
      CHARACTER(LEN=1) :: astring

      ALLOCATE(nm(nlnpe,nlpe,nlev))
      ALLOCATE(np(nlnpe,nlpe,nlev))
      ALLOCATE(cellno(nclprc))

      ALLOCATE(land(nlnpe,nlpe))
      ALLOCATE(o3col(nlnpe,nlpe))
      ALLOCATE(xx(nc,nclprc))
      ALLOCATE(m0(numchem,nlnpe,nlpe,nlev))
      ALLOCATE(conc(nc,nlnpe,nlpe,nlev))
      ALLOCATE(mconc(numchem,nlnpe,nlpe,nlev))
      ALLOCATE(sdconc(numchem,nlnpe,nlpe,nlev))
      ALLOCATE(o3conc(nlnpe,nlpe,nlev))
      ALLOCATE(ch4conc(nlnpe,nlpe,nlev))
      ALLOCATE(estore(nc,nlnpe,nlpe))
      ALLOCATE(dstore(nc,nlnpe,nlpe))
      ALLOCATE(pos(4,nclprc))
      ALLOCATE(o3_ls(nlnpe,nlpe,nmetlev))
      ALLOCATE(so2em(nlnpe,nlpe,nlev))
      ALLOCATE(be7em(nlnpe,nlpe,nlev))
      ALLOCATE(be10em(nlnpe,nlpe,nlev))
      ALLOCATE(lnoxem(nlnpe,nlpe,nlev))
      ALLOCATE(acnoxem(nlnpe,nlpe,nlev))
      ALLOCATE(emiss(nc,nlnpe,nlpe))
      ALLOCATE(nflashes(nlnpe,nlpe,2))

!kk   Get PATH names
      CALL IN_STOCHEM_CHM_GET_PATH

! Initialise grid constants
! DEPENDS ON: ini_lev_const
      CALL INI_LEV_CONST(z_top_of_model,first_constant_r_rho_level)

! Set the number of days in each month.
! DEPENDS ON: set_daym
      CALL SET_DAYM(daym,year)
      ysec = REAL(SUM(daym)*daysec)

      IF (mype==0) THEN
        WRITE(6,*) 'Called Initialise_STOCHEM'
        max_cells=0       ! holds max NFILL from POOL.
      END IF

      IF (mype==0) THEN
        WRITE(6,*) 'Start Year: ',year
        WRITE(6,*) 'Start Month: ',month
        WRITE(6,*) 'Start Day: ',day
      END IF

! Broadcast the date etc. to all PEs
      CALL GC_SSYNC(nproc,info)
      CALL GC_IBCAST(imsg,1,0,nproc,info,year)
      CALL GC_IBCAST(imsg,1,0,nproc,info,month)
      CALL GC_RBCAST(imsg,1,0,nproc,info,day)
      CALL GC_SSYNC(nproc,info)

      WRITE(6,*) 'mype=',mype,'Start : ',day,'/',month,'/',year

! Aircraft ON (NASA data for 1992)  and lightning ON
      IF (scenario == 'pi' .and. lightningon) THEN
        WRITE(6,*)'2. Lightning only NOx'
      ELSE IF (lightningon .and. aircrafton) THEN
        WRITE(6,*)'1. Lightning and Aircraft NOx'
      ELSE IF (lightningon .AND. .NOT. aircrafton) THEN
        WRITE(6,*)'2. Lightning only NOx'
      ELSE IF (.NOT. lightningon .AND. aircrafton) THEN
        WRITE(6,*)'3. Aircraft only NOx'
      ELSE
        WRITE(6,*)'4. No Lightning and Aircraft NOx'
      END IF

!PHOTOLYSIS STARTS AT MIDNIGHT GMT 21st June

! DEPENDS ON: secs
      time = SECS(day,month,year)

! Set up system_seed using standard system values
      CALL FRANVI(system_seed)

! SEED,  For normalised random no. generator;
! SEED2, For selection of random no. from array
      seed2 = 0
      seed3 = 0

! Initialise random number seeds using system_seed
! Shift each seed by mype places.
! Start seed3 halfway through the sequence.
      j = ransize / 2
      DO i = 1, ransize
        seed2(mype,i) = system_seed(1+MOD(i-1+mype,ransize))
        seed3(mype,i) = system_seed(1+MOD(i-1+mype+j,ransize))
      END DO

! Write the number of molecules per air parcel
      WRITE(6,*) 'lmolec: ',lmolec,' molecules'

! Calculate the number of 3D fluxes requested and allocate arrays
      num3dflux = 0
      DO i=1,numflux
        astring = fluxnames(i)(len_flx_str:len_flx_str)
        IF (astring == '1') num3dflux = num3dflux + 1
      END DO

! To avoid memory bank conflicts, make num3dflux_dim odd
      num3dflux_dim = num3dflux
      IF (mod(num3dflux,2) == 0) num3dflux_dim = num3dflux_dim + 1

      WRITE(6,*) num3dflux,' 3D Fluxes Requested'

      ALLOCATE(tflux(num3dflux_dim,nlnpe,nlpe,nlev))
      ALLOCATE(cellflux(num3dflux_dim,nclprc))

      CALL GETLIS

! Initialise air parcel positions
! DEPENDS ON: inipos
      CALL INIPOS(pos,cellno,nfill,orog,z_top_of_model,                 &
     &  first_constant_r_rho_level)

! Calculate number of parcels per 3D Eulerian gridbox
      nnn = 0
! DEPENDS ON: aindex
      CALL AINDEX(nfill,pos,ipos)
      DO j=1,nfill
        nnn(ipos(1,j)-lndat+1,ipos(2,j)-ltdat+1,ipos(3,j)) =            &
     &    nnn(ipos(1,j)-lndat+1,ipos(2,j)-ltdat+1,ipos(3,j))+1
      END DO

! Initialise concentrations of chemical species in each air parcel
! DEPENDS ON: inicon
      CALL INICON(xx,ipos,nfill,month,year)

      tflux = 0.0
      totflu =0.0
      totm0 = 0.0
      totavg =0.0
      daym0 = day
      CALL ZERO_STATS

! Initialise the Eulerian concentrations without mixing
! DEPENDS ON: plot
      CALL PLOT(conc,xx,ipos,nnn,0.,0.,nfill)

! Initialise the ozone profile used for photolysis
      o3conc  = conc(i_o3,:,:,:)
      ch4conc = conc(i_ch4,:,:,:)
      WHERE(o3conc < ppb_vmr_o3) o3conc = o3conc_init
      WHERE(ch4conc < ppb_vmr_ch4) ch4conc = ch4conc_init

! Initialise the land mask
      landum = 0.0
      WHERE(orog > 0.0) landum = 1.0
! DEPENDS ON: met2data
      CALL MET2DATA(land,landum,1,1)

! Read in 72x36xNMETLEV Li & Shine O3 climatology for the month
! used to create o3col
! DEPENDS ON: reado3
      CALL READO3(o3_ls(:,:,1:nmetlev),month,z_top_of_model,            &
     &  first_constant_r_rho_level)

! Calculate a column O3 value (0-model top mb) : o3col (Dobson units)
! DEPENDS ON: o3column
      CALL O3COLUMN(o3_ls,o3col)

! Read in all the monthly emission fields to obtain the annual totals
      Annual_Soil_NOx = 0.0
      Annual_Land_NVOC = 0.0
      Annual_Ocean_NVOC = 0.0
      Annual_Isop = 0.0
      Annual_DMS = 0.0
      DO mm=1,12
        WRITE(smonth,'(I2.2)') mm
        OPEN(22,FILE=TRIM(emdir)//'soilnox'//smonth//'_96_72.dat',      &
     &    STATUS='OLD')
        READ(22,*) temp_emiss
        CLOSE(22)
        Annual_Soil_NOx=Annual_Soil_NOx+SUM(temp_emiss)
        OPEN(22,FILE=TRIM(emdir)//'nvoc_land'//smonth//'_96_72.dat',    &
     &    STATUS='OLD')
        READ(22,*) temp_emiss
        CLOSE(22)
        Annual_Land_NVOC=Annual_Land_NVOC+SUM(temp_emiss)
        OPEN(22,FILE=TRIM(emdir)//'nvoc_ocean'//smonth//'_96_72.dat',   &
     &    STATUS='OLD')
        READ(22,*) temp_emiss
        CLOSE(22)
        Annual_Ocean_NVOC=Annual_Ocean_NVOC+SUM(temp_emiss)
        OPEN(22,FILE=TRIM(emdir)//'isop'//smonth//'_96_72.dat',         &
     &    STATUS='OLD')
        READ(22,*) temp_emiss
        CLOSE(22)
        Annual_Isop=Annual_Isop+SUM(temp_emiss)
        OPEN(22,FILE=TRIM(emdir)//'dms_ocean'//smonth//'_96_72.dat',    &
     &    STATUS='OLD')
        READ(22,*) temp_emiss
        CLOSE(22)
        Annual_DMS=Annual_DMS+SUM(temp_emiss)
      END DO

      annual_soil_nox = annual_soil_nox / 12.0
      annual_land_nvoc = annual_land_nvoc / 12.0
      annual_ocean_nvoc = annual_ocean_nvoc / 12.0
      annual_isop = annual_isop / 12.0
      annual_dms = annual_dms / 12.0

      WRITE(6,'(a30,e13.5,a8)') 'Annual_Soil_NOx Emissions: ',          &
     &  Annual_Soil_NOx,' (kg/yr)'
      WRITE(6,'(a30,e13.5,a8)') 'Annual_Land_NVOC Emissions: ',         &
     &  Annual_Land_NVOC,' (kg/yr)'
      WRITE(6,'(a30,e13.5,a8)') 'Annual_Ocean_NVOC Emissions: ',        &
     &  Annual_Ocean_NVOC,' (kg/yr)'
      WRITE(6,'(a30,e13.5,a8)') 'Annual_Isop Emissions: ',              &
     &  Annual_Isop,' (kg/yr)'
      WRITE(6,'(a30,e13.5,a8)') 'Annual_DMS Emissions: ',               &
     &  Annual_DMS,' (kg/yr)'

      emiss = 0.0
      estore =0.0
      ddepo = 0.0
      dstore =0.0

! Calculate emissions
!   EMISS has units (molecules s^-1 per grid square)
!                   (------------------------------)
!                   (     molecules per cell       )

! Call appropriate emissions routine.
      SELECT CASE(scenario)
      CASE('fi')
! DEPENDS ON: emcalc_a1fi
        CALL EMCALC_A1FI(emiss,month,year)
      CASE('ab')
! DEPENDS ON: emcalc_a1b
        CALL EMCALC_A1B(emiss,month,year)
      CASE('a2')
! DEPENDS ON: emcalc_a2
        CALL EMCALC_A2(emiss,month,year)
      CASE('b1')
! DEPENDS ON: emcalc_b1
        CALL EMCALC_B1(emiss,month,year)
      CASE('b2')
! DEPENDS ON: emcalc_b2
        CALL EMCALC_B2(emiss,month,year)
      CASE('bu')
! DEPENDS ON: emcalc_cle
        CALL EMCALC_CLE(emiss,month,year)
      CASE('mf')
! DEPENDS ON: emcalc_mfr
        CALL EMCALC_MFR(emiss,month,year)
      CASE('A2')
! DEPENDS ON: emcalc_a2_
        CALL EMCALC_A2_(emiss,month,year)
      CASE('pi')
! DEPENDS ON: emcalc_pi
        CALL EMCALC_PI(emiss,month,year)
      CASE DEFAULT
        cmessage = 'No emission scenario selected'
        WRITE(6,*) cmessage
! DEPENDS ON: ereport
        CALL EREPORT('Initialise_Stochem',1,cmessage)
      END SELECT

! Read in and convert NO2 emissions from aircraft
! DEPENDS ON: lightread
      CALL LIGHTREAD(acnoxem,so2em,be7em,be10em,month,year,             &
     &  totacnoxem,totso2em,totbe7em,totbe10em,p,lnp)

! Calculate Lightning NOx emissions interactively
      IF (lightningon) THEN
! DEPENDS ON: lightnox
        CALL LIGHTNOX(cc_base_level,cc_top_level,t,acr,orog,            &
     &    land_fraction,totlnoxem,lnoxem,nflashes,z_top_of_model,       &
     &    first_constant_r_rho_level)
      ELSE
        lnoxem = 0.0
        totlnoxem = 0.0
        nflashes = 0.0
      END IF

! Calculate isoprene and monoterpene emissions
! DEPENDS ON: emcalc_natural_voc
      CALL EMCALC_NATURAL_VOC(month,day,totpar,dirpar,gsf,t0tile,       &
     &  lai_ft,emiss(i_c5h8,:,:),terpem)

      CALL RESTART               ! look for a STOCHEM dump file

      StochDumpRequest = .false. ! no dump at start

      END SUBROUTINE INITIALISE_STOCHEM
!#######################################################################

!#######################################################################
      SUBROUTINE MAKE_SWAP_MAP
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Even out cell distribution between processors
!-                         Internal procedure to STOCHEM
!-
!-   Inputs  :
!-   Outputs : CELLSWAP - indirectly
!-   Controls:
!-
! Current Code Owner: M.G. Sanderson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   08/04/05   Created. R. Johanni
!-
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------

      INTEGER :: naverage, nremain, num_send, num_recv, npe_send, n

! Calculate integer <= mean no. of cells per PE
      naverage=ncell/nproc
! But NPROC doesn't necessarily divide NCELL exactly.
      nremain=ncell-naverage*nproc   ! calculate remainder
      nc_actual=0
      nc_actual(mype)=nfill ! number of cells on this PE

! Sum nc_actual over all PEs, equivalent to broadcasting all
! non zero elements.
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(nproc,info)
      CALL GC_ISUM(nproc,nproc,info,nc_actual)
      CALL GC_SSYNC(nproc,info)

      IF (mype == 0) THEN
        DO n=0,nproc-1
          max_cells(n)=MAX(max_cells(n),nc_actual(n))
        END DO
      END IF

! Calculate how many cells each PE should have

      nc_target(:) = naverage
! the first nremain PEs get an extra cell
      DO i=0,nremain-1
        nc_target(i) = nc_target(i)+1
      END DO

      WRITE(6,*) 'MAKE_SWAP_MAP nc_actual: ',nc_actual(0:nproc-1)
      WRITE(6,*) 'MAKE_SWAP_MAP nc_target: ',nc_target(0:nproc-1)

! For safety only, should never happen

      IF(SUM(nc_actual(0:nproc-1)) /= SUM(nc_target(0:nproc-1))) THEN
        cmessage=' Error pooling cells: nc_actual != nc_target'
        write(6,*) cmessage
        write(6,*) 'SUM(nc_actual):',SUM(nc_actual(0:nproc-1))
        write(6,*) 'SUM(nc_target):',SUM(nc_target(0:nproc-1))
! DEPENDS ON: ereport
        CALL EREPORT('POOL',1,cmessage)
        CALL abort
      ENDIF

! Now produce a map which PEs get how many cells from whom
! we do this by assigning one or more senders to every receiver

      npe_send = -1 ! Next PE which sends cells
      num_send = 0  ! number of cells remaining on npe_send
      nmap = 0      ! number of map entries

      DO n=0,nproc-1

! check if PE n is a receiver:

        IF (nc_actual(n) >= nc_target(n)) CYCLE
        num_recv = nc_target(n)-nc_actual(n)

! search senders

        DO WHILE(num_recv > 0)

          IF (num_send==0) THEN
! search the next sender
            DO WHILE (npe_send<nproc-1)
              npe_send = npe_send + 1
              IF (nc_actual(npe_send) > nc_target(npe_send)) THEN
! found one
                num_send = nc_actual(npe_send) - nc_target(npe_send)
                EXIT
              END IF
            END DO
! For safety only, should never happen:
            IF (num_send == 0) THEN
! left the above loop without finding a sender
              cmessage=' Error pooling cells: no sender'
              WRITE(6,*) cmessage
! DEPENDS ON: ereport
              CALL EREPORT('POOL',1,cmessage)
              CALL abort
            END IF
          END IF

! add map entry

          nmap = nmap+1
          balance_map(1,nmap) = n        ! Receiver
          balance_map(2,nmap) = npe_send ! Sender

          IF(num_send >= num_recv) THEN
            balance_map(3,nmap) = num_recv ! Number of cells
            num_send = num_send - num_recv
            num_recv = 0
          ELSE
            balance_map(3,nmap) = num_send ! Number of cells
            num_recv = num_recv - num_send
            num_send = 0
          END IF
        END DO

      END DO

! For safety only, should never happen

      IF (num_send/=0 .OR. num_recv/=0) THEN
        cmessage=' Error pooling cells: cells remain'
        WRITE(6,*) cmessage
! DEPENDS ON: ereport
        CALL EREPORT('POOL',1,cmessage)
        CALL abort
      END IF

      END SUBROUTINE MAKE_SWAP_MAP
      SUBROUTINE MAP1(SEND_MAP,RECV_MAP,FROM,TO,NEWCELL,NDIM,           &
     &  N_ITEMS_SEND,N_ITEMS_RECV)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Set up send and receive maps.
!-                         Internal procedure to STOCHEM
!-
!-   Inputs  : FROM,TO,NDIM,N_ITEMS_SEND,N_ITEMS_RECV
!-   Outputs : SEND_MAP,RECV_MAP
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.4    20/11/97  Created.  W.J. Collins
!  4.5    01/03/98  Reduced size of cell arrays.  W.J. Collins
!  6.1    20/10/04  No change.
!
!-
!VVV  V2.2  MAP1 20/X/99
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, DIMENSION(NCLPRC), INTENT(IN) :: FROM,TO
      INTEGER, DIMENSION(NCLPRC*NPROC), INTENT(IN) :: NEWCELL
      INTEGER, INTENT(IN) :: NDIM,N_ITEMS_SEND,N_ITEMS_RECV
      INTEGER, DIMENSION(S_STRIDE_IN_RECV_ARRAY,NCLPRC), INTENT(INOUT)  &
     &   :: SEND_MAP
      INTEGER, DIMENSION(R_STRIDE_IN_SEND_ARRAY,NCLPRC), INTENT(INOUT)  &
     &   :: RECV_MAP

      SEND_MAP(S_BASE_ADDRESS_IN_SEND_ARRAY,1:N_ITEMS_SEND)=            &
     &  (FROM(1:N_ITEMS_SEND)-CELLBASE-1)*NDIM+1
      SEND_MAP(S_NUMBER_OF_ELEMENTS_IN_ITEM,1:N_ITEMS_SEND)=1
      SEND_MAP(S_STRIDE_IN_SEND_ARRAY,      1:N_ITEMS_SEND)=NDIM
      SEND_MAP(S_ELEMENT_LENGTH,            1:N_ITEMS_SEND)=NDIM
      SEND_MAP(S_BASE_ADDRESS_IN_RECV_ARRAY,1:N_ITEMS_SEND)=            &
     &  (NEWCELL(FROM(1:N_ITEMS_SEND))-1)*NDIM+1
      SEND_MAP(S_STRIDE_IN_RECV_ARRAY,      1:N_ITEMS_SEND)=NDIM

      RECV_MAP(R_BASE_ADDRESS_IN_RECV_ARRAY,1:N_ITEMS_RECV)=            &
     &  (NEWCELL(TO(1:N_ITEMS_RECV))-1)*NDIM+1
      RECV_MAP(R_NUMBER_OF_ELEMENTS_IN_ITEM,1:N_ITEMS_RECV)=1
      RECV_MAP(R_STRIDE_IN_RECV_ARRAY,      1:N_ITEMS_RECV)=NDIM
      RECV_MAP(R_ELEMENT_LENGTH,            1:N_ITEMS_RECV)=NDIM
      RECV_MAP(R_BASE_ADDRESS_IN_SEND_ARRAY,1:N_ITEMS_RECV)=            &
     &  MOD(TO(1:N_ITEMS_RECV)-1,NCLPRC)*NDIM+1
      RECV_MAP(R_STRIDE_IN_SEND_ARRAY,      1:N_ITEMS_RECV)=NDIM
      SEND_MAP(:,N_ITEMS_SEND+1:)=0
      RECV_MAP(:,N_ITEMS_RECV+1:)=0

      END SUBROUTINE MAP1
! ######################################################################
      SUBROUTINE MONTHLY_OUTPUT
!-----------------------------------------------------------------------
!
!     Produce output as month changes, archive the files
!     Internal procedure to STOCHEM.
!
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.1    14/03/00  Created.  D.S. Stevenson
!  5.1    05/09/00  Archives station files.  C.E. Johnson
!  5.3    25/01/01  daym0, daymas added.  C.E. Johnson
!  6.1    21/10/04  No change.
!
!
!VVV  V2.4  monthly_output 25/I/00  added daym0,daymas
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      REAL :: daymas
!
! Output current set of statistics,
! DEPENDS ON: mcalc
      CALL MCALC(mass,totmas,xx,clist,ipos,pos,tropz,nchem,nfill)
      daymas=day
! DEPENDS ON: calsta
      CALL CALSTA(sdconc,mconc,np,nph,clist,nchem)

! Sum TOTs over all processors
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(nproc,info)
      CALL GC_RSUM(2*nc,nproc,info,totavg)
      CALL GC_RSUM(2*nc,nproc,info,totm0)
      CALL GC_RSUM(2*nc,nproc,info,totmas)
      CALL GC_RSUM(numflux,nproc,info,totflu)
      info=GC_SHM_PUT
      CALL GC_SSYNC(nproc,info)

! Output global inventories and fluxes, 3-D fluxes and
! the 3-D species mixing ratios.
! DEPENDS ON: outbud
      CALL OUTBUD(nm,mass,m0,tflux,clist,nchem,fnames,nflux,            &
     &  totm0,daym0,totmas,daymas,totavg/navg,totflu,navg,flist,        &
     &  stash_code,month,year,mconc,fixhd12,period,umstepno,            &
     &  z_top_of_model,first_constant_r_rho_level)

      IF (mype == 0) THEN
! Archive the monthly mean pp file.
! DEPENDS ON: archive
        CALL ARCHIVE('p',filetype2,month,year)
! Archive the monthly inventory and budget file.
! DEPENDS ON: archive
        CALL ARCHIVE('t',filetype2,month,year)
! Archive the monthly station data file
! DEPENDS ON: archive
        CALL ARCHIVE('s','m',month,year)
! Archive the cell following data file
! DEPENDS ON: archive
        CALL ARCHIVE('f','m',Month,year)
      ENDIF

      END SUBROUTINE MONTHLY_OUTPUT
!#######################################################################
#ifdef OLDSWAP2
      SUBROUTINE POOL
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Even out cell distribution between processors
!-                         Internal procedure to STOCHEM
!-
!-   Inputs  :
!-   Outputs : CELLSWAP - indirectly
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.4    30/09/97  Created.  W.J. Collins
!  4.5    01/03/98  Reduced size of cell arrays.  W.J. Collins
!  6.1    21/10/04  No change.
!
!-
!VVV  V2.2  POOL 20/X/99
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, DIMENSION(0:maxproc-1) :: numcell
      INTEGER :: i,j,k,naverage,nremain,nprev,nextra,                   &
     &  nadd,extra,swapsize
      INTEGER, DIMENSION(mxcoll) :: clswp1
      COMMON /CPOOL/numcell,clswp1 ! aligns variables for shared memory

! Calculate integer <= mean no. of cells per PE
      naverage=ncell/nproc
! But NPROC doesn't necessarily divide NCELL exactly.
      nremain=ncell-naverage*nproc   ! calculate remainder
      cellswap=0
      numcell=0
      numcell(mype)=nfill ! number of cells on this PE

! Sum NUMCELL over all PEs, equivalent to broadcasting all
! non zero elements.
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(nproc,info)
      CALL GC_ISUM(maxproc,nproc,info,numcell)
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(nproc,info)

      IF (mype == 0) THEN
        DO i=0,nproc-1
          max_cells(i)=MAX(max_cells(i),numcell(i))
        END DO
      END IF

! NREMAIN processors are allowed an extra cell to make up the total
      extra=0
! Count number of processors 0:mype-1 that have already got an extra cel
      IF(COUNT(numcell(0:mype-1)>=naverage+1)<nremain) extra=1

! If too many cells on this processor
      IF(nfill>naverage+extra) THEN
        i=0
        DO j=1,nfill
! flag cells to be swapped
          cellswap(1,j+cellbase)=mype+1
          i=i+1
! Have we swapped enough cells?
          IF(i>=nfill-naverage-extra) EXIT
        END DO
        IF(j>nfill) THEN
          cmessage=' Error pooling cells: J > NFILL'
          write(6,*) cmessage
! DEPENDS ON: ereport
          CALL EREPORT('POOL',1,cmessage)
        ENDIF
      END IF

      swapsize=nclprc*nproc
! Sum CELLSWAP(1) over all PEs, equiv. to broadcasting all
! non zero elements.
      info=GC_SHM_PUT
      DO k=0,swapsize/mxcoll
        clswp1=0
        clswp1(1:MIN(swapsize-k*mxcoll,mxcoll))=                        &
     &    cellswap(1,k*mxcoll+1:min((k+1)*mxcoll,swapsize))
        CALL GC_SSYNC(nproc,info)
        CALL GC_ISUM(mxcoll,nproc,info,clswp1)
        CALL GC_SSYNC(nproc,info)
        cellswap(1,k*mxcoll+1:MIN((k+1)*mxcoll,swapsize))=              &
     &    clswp1(1:MIN(swapsize-k*mxcoll,mxcoll))
      END DO
!      PRINT *,'NO. CELLS SWAPPED=',COUNT(CELLSWAP(1,:)>0)

! CELLSWAP(1) is now set for all cells to be swapped

      extra=0
! If not enough processors have one extra cell, then we need
! to swap some in.
! no. initially having 1 extra
      nextra=COUNT(numcell(:)  >=naverage+1)
! no. extras added to PEs 0:mype-1
      nadd=COUNT(numcell(0:mype-1)<naverage+1)
      IF(nextra+nadd<nremain) EXTRA=1

! Too few cells, need to add some.
      IF(numcell(mype)<naverage+extra) THEN
! number of cells added to PEs 0:mype-1, allowing for extras
        nprev=SUM(MAX(naverage-numcell(0:mype-1),0))
! (no. of extras)
        IF(nextra<nremain) nprev=nprev+MIN(nremain-nextra,nadd)
        i=0
        DO j=1,nclprc*nproc
          IF(CELLSWAP(1,j)/=0) THEN ! found a cell to swap
            i=i+1
! Check whether the cell will have been added to a previous PE
            IF(i>nprev) cellswap(2,j)=mype+1
! added enough
            IF(i>=nprev+naverage+extra-numcell(mype)) EXIT
          END IF
        END DO
        IF(j>nclprc*nproc) THEN
          cmessage=' Error assigning cells'
! DEPENDS ON: ereport
          CALL EREPORT('POOL',1,cmessage)
        END IF
      END IF
! Zero CELLSWAP if not swapping to this PE
      WHERE(cellswap(2,:)==0) cellswap(1,:)=0

! Don't swap a cell to the same processor
      WHERE(cellswap(1,:)==cellswap(2,:).and.cellswap(1,:)>0)
        cellswap(1,:)=0
        cellswap(2,:)=0
      END WHERE

      END SUBROUTINE POOL
#endif
!#######################################################################
      SUBROUTINE RESTART
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Restarts from dumped variables. Internal
!-                         procedure to STOCHEM
!-
!-   Inputs  :
!-   Outputs : Indirectly- XX,POS,ESTORE,DSTORE,TFLUX,TOTFLU,FLIST,
!-             FNAMES,CONC,MCONC,SDCONC,O3CONC,NM,NP,
!-             M0,TIME,YEAR,MONTH,DAY,SEED2,SEED3,NFLUX
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.4    17/01/95  Created.  W.J. Collins
!  4.5    17/01/96  Added TOTAVG, NAVG.  W.J. Collins
!  4.5    12/06/98  Moved time etc. to start of dump.  W.J. Collins
!  4.5    12/10/98  Added O3CONC to dump.  D.S. Stevenson
!  4.5    20/04/99  Stop coupled run if start attempted mid run.
!                   C.E. Johnson
!  5.0    11/05/99  UM-style filenames. C.E. Johnson
!  5.0    17/01/96  Added DAYM0. C.E. Johnson
!  5.5    10/10/03  Uses new version of AINDEX. Seed array zeroing bug
!                   fixed. Added re_num3dflux_dim. M.G. Sanderson
!  6.1    23/09/04  Reads in integer random number seeds. M.G. Sanderson
!

!VVV  V5.2  RESTART  2/8/01 - B deleted
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER :: ierr
      INTEGER :: info
      INTEGER :: imsg=1
      INTEGER :: i
      INTEGER :: j
      INTEGER :: iproc
      INTEGER :: icell
      INTEGER :: icellno
      INTEGER :: re_month
      INTEGER :: re_year
      INTEGER :: re_nproc
      INTEGER :: re_nlong
      INTEGER :: re_mnlat
      INTEGER :: re_nlev
      INTEGER :: re_ncell
      INTEGER :: re_nc
      INTEGER :: re_numflux
      INTEGER :: re_num3dflux
      INTEGER :: re_num3dflux_dim
      INTEGER :: re_numchem
      INTEGER :: icode
      INTEGER :: irestart

      REAL    :: re_day
      REAL    :: re_time
      REAL, DIMENSION(4)  :: rpos
      REAL, DIMENSION(NC) :: rxx

      CHARACTER(LEN=8), DIMENSION(numchem) :: re_outnames
      CHARACTER(LEN=len_flx_str), DIMENSION(numflux) :: re_fluxnames
      CHARACTER(LEN=72) :: cmessage
      CHARACTER(LEN=14) :: filename
      COMMON /REST/ ierr

      icode = 0
      irestart = 0

! Read in cell data on all PEs. File now has current time.
! DEPENDS ON: findname
      CALL FINDNAME('d','a','c',0,0,filename)
      OPEN(52,FILE=filename,STATUS='OLD',                               &
     &  FORM='UNFORMATTED',IOSTAT=IERR)

! Test to see if file exists.
      IF (ierr == 0) THEN
        IF (mype == 0) WRITE(6,*) 'DUMP FILE FOUND,',filename
! Read in dump common info and test for compatibility with time
! and chemical model dimensions etc.
        nfill = 0
        READ(52) re_time,re_year,re_month,re_day
        IF (mype==0) WRITE(6,*) 'Reading dump file from: ',re_day,'-',  &
     &    re_month,'-',re_year
        IF (re_time /= time) THEN
          icode = -1
          cmessage = 'STOCHEM dump time incorrect'
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
          WRITE(6,*) cmessage,' time,re_time: ',time,re_time
          WRITE(6,*) 'STOCHEM dump has date:',                          &
     &      ' DAY: ',re_day,' MONTH: ',re_month,' YEAR: ',re_year
        END IF

        READ(52) re_nproc,re_nlong,re_mnlat,re_nlev,re_ncell,re_nc
        IF (re_nproc /= nproc) THEN
          icode = 1
          cmessage='STOCHEM dump incompatible with number of PEs'
          WRITE(6,*) cmessage,' RE_NPROC: ',re_nproc
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',ICODE,cmessage)
        END IF
        IF (re_nlong /= nlong) THEN
          icode = 1
          cmessage='STOCHEM dump has incompatible NLONG'
          WRITE(6,*) cmessage,' RE_NLONG: ',re_nlong
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
        END IF
        IF (re_mnlat /= mnlat) THEN
          icode = 1
          cmessage = 'STOCHEM dump has incompatible MNLAT'
          WRITE(6,*) cmessage,' RE_MNLAT: ',re_mnlat
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
        END IF
        IF (re_nlev /= nlev) THEN
          icode = 1
          cmessage = 'STOCHEM dump has incompatible NLEV'
          WRITE(6,*) cmessage,' RE_NLEV: ',re_nlev
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
        END IF
        IF (re_ncell /= ncell) THEN
          icode = 1
          cmessage = 'STOCHEM dump has incompatible NCELL'
          WRITE(6,*) cmessage,' RE_NCELL: ',re_ncell
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
        END IF
        IF (re_nc /= nc) THEN
          icode = 1
          cmessage = 'STOCHEM dump has incompatible NC'
          WRITE(6,*) cmessage,' RE_NC: ',re_nc
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
        END IF

        READ(52) re_numflux,re_num3dflux,re_num3dflux_dim,              &
     &    re_numchem,re_outnames,re_fluxnames
        IF (re_numflux /= numflux) THEN
          icode = -1
          irestart = 1
          cmessage = 'STOCHEM dump has incompatible NUMFLUX'
          WRITE(6,*) cmessage,' RE_NUMFLUX: ',re_numflux
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
        END IF
        IF (re_num3dflux /= num3dflux) THEN
          icode = -1
          irestart = 1
          cmessage = 'STOCHEM dump has incompatible NUM3DFLUX'
          WRITE(6,*) cmessage,' RE_NUM3DFLUX: ',re_num3dflux
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
        END IF
        IF (re_num3dflux_dim /= num3dflux_dim) THEN
          icode = -1
          irestart = 1
          cmessage = 'STOCHEM dump has incompatible NUM3DFLUX_DIM'
          WRITE(6,*) cmessage,' RE_NUM3DFLUX_DIM: ',re_num3dflux_dim
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
        ENDIF
        IF (re_numchem /= numchem) THEN
          icode = -1
          irestart = 1
          cmessage = 'STOCHEM dump has incompatible NUMCHEM'
          WRITE(6,*) cmessage,' RE_NUMCHEM: ',re_numchem
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
        END IF

        IF (irestart == 0) THEN
! OUTNAMES and FLUXNAMES have the correct lengths, but we need to
! check that output species and fluxes are the same ones, and
! set IRESTART accordingly.
          DO i=1,numchem
            IF (re_outnames(i) /= outnames(i)) THEN
              icode = -1
              irestart = 1
              cmessage = 'STOCHEM dump has incompatible OUTNAMES'
              WRITE(6,*) cmessage,' RE_OUTNAMES: ',re_outnames
! DEPENDS ON: ereport
              CALL EREPORT('RESTART',icode,cmessage)
            END IF
          END DO
          DO i=1,numflux
            IF (re_fluxnames(i) /= fluxnames(i)) THEN
              icode = -1
              irestart = 1
              cmessage = 'STOCHEM dump has incompatible FLUXNAMES'
              WRITE(6,*) cmessage,' RE_FLUXNAMES: ',re_fluxnames
! DEPENDS ON: ereport
              CALL EREPORT('RESTART',ICODE,cmessage)
            END IF
          END DO
        END IF

        IF (irestart /= 0) THEN
! Fluxes or output species are incompatible, but restart anyway
          icode = -1
          cmessage='STOCHEM restarting, but some output info missing'
          WRITE(6,*) cmessage
! DEPENDS ON: ereport
          CALL EREPORT('RESTART',icode,cmessage)
          CLOSE(52)
! Do a restricted restart - some output info is missing, so
! the next months output will be incomplete.
!         CALL RESTART_R(re_numflux,re_numchem,re_num3dflux)
        ELSE
! If all is well, read rest of dump
          DO j=1,ncell
            READ(52) iproc,icell,rpos,rxx,icellno
            IF (iproc == mype) THEN
              nfill = nfill + 1
              pos(:,nfill) = rpos
              xx(:,nfill) = rxx
              cellno(nfill) = icellno
            END IF
          END DO
! Read in Eulerian data on all PEs then set smaller arrays as appropriat
! Slightly extravagant
          DO i=1,nc
            READ(52) global(:,:,1)
            estore(i,:,:)=global(lndat:lndat+nlnpe-1,                   &
     &                           ltdat:ltdat+nlpe-1,1)
            READ(52) global(:,:,1)
            dstore(i,:,:)=global(lndat:lndat+nlnpe-1,                   &
     &                           ltdat:ltdat+nlpe-1,1)
            READ(52) global
            conc(i,:,:,:)=global(lndat:lndat+nlnpe-1,                   &
     &                           ltdat:ltdat+nlpe-1,:)
          END DO
          DO i=1,numchem
            READ(52) global
            mconc(i,:,:,:)=global(lndat:lndat+nlnpe-1,                  &
     &                            ltdat:ltdat+nlpe-1,:)
            READ(52) global
            sdconc(i,:,:,:)=global(lndat:lndat+nlnpe-1,                 &
     &                             ltdat:ltdat+nlpe-1,:)
            READ(52) global
            m0(i,:,:,:)=global(lndat:lndat+nlnpe-1,                     &
     &                         ltdat:ltdat+nlpe-1,:)
          END DO
          DO i=1,num3dflux_dim
            READ(52) global
            tflux(i,:,:,:)=global(lndat:lndat+nlnpe-1,                  &
     &                            ltdat:ltdat+nlpe-1,:)
          END DO
          READ(52) global
          o3conc=global(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1,:)
          READ(52) global
          nm=INT(global(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1,:))
          READ(52) global
          np=INT(global(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1,:))
! Read in all other data on PE 0
          IF (mype==0) THEN
            READ(52) totflu
            READ(52) flist,fnames
            READ(52) stash_code
            READ(52) seed2
            READ(52) seed3
            READ(52) nflux,nstat
            READ(52) totm0
            READ(52) totavg
            READ(52) navg
            READ(52) daym0
            READ(52) annual_soil_nox
            READ(52) annual_land_nvoc
            READ(52) annual_ocean_nvoc
            READ(52) annual_isop
            READ(52) annual_dms
            WRITE(6,*) 'RESTART FILE EXISTS, VARIABLES RESTORED'
            WRITE(6,*) 'TIME: ',time
            WRITE(6,*) 'DAY: ',day,'MONTH: ',month,'YEAR: ',year
          END IF
          CLOSE(52)
        END IF      ! IRESTART
      ELSE
        IF (mype==0) WRITE(6,*) 'RESTART DUMP FILE NOT FOUND,',         &
     &   filename,' STARTING FROM SCRATCH'
! Inhibit the model from starting from scratch unless at start of UM run
         IF (UMstepNo > 20) then
           cmessage = 'STOCHEM restart dump file not found.'
           WRITE(6,*) cmessage,'FILENAME: ',filename
! DEPENDS ON: ereport
           CALL EREPORT('RESTART',1,cmessage)
         END IF
      END IF

      IF (ierr == 0) THEN ! if dump file exists
        WRITE(6,*) 'RESTART FILE EXISTS, VARIABLES RESTORED'
! Broadcast common infomation from PE 0 to all other PEs
        CALL GC_SSYNC(nproc,info)
        CALL GC_IBCAST(imsg,1,0,nproc,info,nflux)
        CALL GC_IBCAST(imsg,1,0,nproc,info,nstat)
        CALL GC_IBCAST(imsg,1,0,nproc,info,navg)
        CALL GC_IBCAST(imsg,2*numflux,0,nproc,info,flist)
        CALL GC_IBCAST(imsg,numflux,0,nproc,info,stash_code)
        CALL GC_IBCAST(imsg,maxproc*ransize,0,nproc,info,seed2)
        CALL GC_IBCAST(imsg,maxproc*ransize,0,nproc,info,seed3)
        INFO = GC_SHM_PUT
        CALL GC_SSYNC(nproc,info)

! Set the values of SEEDs to zero that aren't for this PE.
        DO i = 0, maxproc-1
          IF (i /= mype) THEN
            seed2(i,:) = 0
            seed3(i,:) = 0
          END IF
        END DO

! Set up NNN

        nnn=0
! DEPENDS ON: aindex
        CALL AINDEX(nfill,pos,ipos)
        DO j=1,nfill
          nnn(ipos(1,j)-lndat+1,ipos(2,j)-ltdat+1,ipos(3,j))=           &
     &      nnn(ipos(1,j)-lndat+1,ipos(2,j)-ltdat+1,ipos(3,j))+1
        END DO
      END IF

      END SUBROUTINE RESTART
#if defined(OLDSWAP1) || defined(OLDSWAP2)
      SUBROUTINE SWAP
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Swaps cell information between processors.
!-                         Internal procedure to STOCHEM.
!-
!-   Inputs  : also CELLSWAP, POS and XX
!-   Outputs : also POS and XX
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.4    25/09/97  Created.  W.J. Collins
!  4.5    01/03/98  Reduced size of cell arrays. W.J. Collins
!  5.5    29/01/04  Cellflux uses num3dflux_dim dimension. M.G.Sanderson
!  5.5    11/03/04  Now calls new routine kk_ralltoalle instead of
!                   gcg_ralltoalle as former is faster. K. Ketelsen
!  6.1    21/10/04  No change.
!-
!VVV  V2.2  SWAP 20/X/99
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, DIMENSION(nclprc) :: hole
      INTEGER            :: sarr_len
      INTEGER            :: rarr_len
      INTEGER, PARAMETER :: unused = 0
      INTEGER, PARAMETER :: gid = unused
      INTEGER, PARAMETER :: flag =unused
      INTEGER :: info
      INTEGER :: n_items_send
      INTEGER :: n_items_recv
      INTEGER :: ndim
      INTEGER, DIMENSION(s_stride_in_recv_array, nclprc) :: send_map
      INTEGER, DIMENSION(r_stride_in_send_array, nclprc) :: recv_map
      INTEGER, DIMENSION(nclprc) :: from,to
      INTEGER, DIMENSION(nclprc*nproc) :: newcell
      integer :: j
      integer :: jj
      integer :: k
      integer :: swapsize
      REAL, DIMENSION(nclprc) :: rcellno
      INTEGER, DIMENSION(mxcoll) :: clswp1
      INTEGER, DIMENSION(mxcoll) :: nwswp1
!      COMMON /CSWAP/ CLSWP1,NWSWP1,RCELLNO

      n_items_send=0
      n_items_recv=0
      newcell=0
      swapsize=nclprc*nproc
      rcellno=REAL(cellno)

! Sum CELLSWAP over all PEs, equivalent to broadcasting all
! non zero elements.
!      Done in two sections to avoid overflow of workspace.
      DO j=1,2
        info=gc_shm_put
        DO k=0,swapsize/mxcoll
          clswp1=0
          clswp1(1:min(swapsize-k*mxcoll,mxcoll))=                      &
     &      cellswap(j,k*mxcoll+1:min((k+1)*mxcoll,swapsize))
          CALL gc_ssync(nproc,info)
          CALL gc_isum(mxcoll,nproc,info,clswp1)
          CALL gc_ssync(nproc,info)
          cellswap(j,k*mxcoll+1:min((k+1)*mxcoll,swapsize))=            &
     &      clswp1(1:min(swapsize-k*mxcoll,mxcoll))
        END DO
      END DO

      HOLE=0
! Have to remember to subtract 1 from CELLSWAP each time to get
! correct PE.
      DO j=1,nclprc*nproc
        IF(cellswap(1,j)-1==mype) THEN
          n_items_send=n_items_send+1
          send_map(s_destination_pe,n_items_send)=cellswap(2,j)-1
          from(n_items_send)=j ! from = list of cells sent from this pe
          hole(n_items_send)=j-cellbase ! cell no longer on this pe
        END IF

        IF(cellswap(2,j)-1==mype) THEN
          n_items_recv=n_items_recv+1
          recv_map(r_source_pe,n_items_recv)=cellswap(1,j)-1
          to(n_items_recv)=j ! TO is list of cells sent to this PE
          newcell(j)=nfill+n_items_recv ! New no. for cell on this PE
        END IF
      ENDDO

! Check if there is enough space allocated
      IF(nfill+n_items_recv > nclprc) THEN
        cmessage='STOCHEM Routine SWAP allocated too many cells '
        WRITE(6,*) cmessage
        WRITE(6,*) 'NFILL: ',NFILL,' N_ITEMS_RECV: ',N_ITEMS_RECV
        WRITE(6,*) 'NCLPRC: ',NCLPRC
! DEPENDS ON: ereport
        CALL EREPORT('SWAP',1,cmessage)
      ENDIF


! Sum NEWCELL over all PEs, equivalent to broadcasting all
! non zero elements.
      info=gc_shm_put
      DO k=0,swapsize/mxcoll
        nwswp1=0
        nwswp1(1:min(swapsize-k*mxcoll,mxcoll))=                        &
     &    newcell(k*mxcoll+1:min((k+1)*mxcoll,swapsize))
        CALL gc_ssync(nproc,info)
        CALL gc_isum(mxcoll,nproc,info,nwswp1)
        CALL gc_ssync(nproc,info)
        newcell(k*mxcoll+1:min((k+1)*mxcoll,swapsize))=                 &
     &    nwswp1(1:min(swapsize-k*mxcoll,mxcoll))
      END DO
! Get send and recv maps. Send using RALLTOALLE

      ndim=4 ! 1st dimension of POS
      CALL map1(send_map,recv_map,from,to,newcell,ndim,                 &
     &  n_items_send,n_items_recv)
      sarr_len=nclprc*ndim
      rarr_len=nclprc*ndim

!  send out all the swapped elements of POS
      info=gc_shm_put
      CALL gc_ssync(nproc,info)
      CALL kk_ralltoalle(pos,send_map,n_items_send,sarr_len,            &
     &                    pos,recv_map,n_items_recv,rarr_len,           &
     &                    gid,flag,info)
      CALL gc_ssync(nproc,info)

      ndim=nc ! 1st dimension of XX
      CALL map1(send_map,recv_map,from,to,newcell,ndim,                 &
     &  n_items_send,n_items_recv)
      sarr_len=nclprc*ndim
      rarr_len=nclprc*ndim

! send out all the swapped elements of XX
      info=gc_shm_put
      CALL gc_ssync(nproc,info)
      CALL kk_ralltoalle(xx,send_map,n_items_send,sarr_len,             &
     &                    xx,recv_map,n_items_recv,rarr_len,            &
     &                    gid,flag,info)
      CALL gc_ssync(nproc,info)

      ndim=num3dflux_dim ! 1st dimension of CELLFLUX
      CALL map1(send_map,recv_map,from,to,newcell,ndim,                 &
     &  n_items_send,n_items_recv)
      sarr_len=nclprc*ndim
      rarr_len=nclprc*ndim

! send out all the swapped elements of CELLFLUX
      info=gc_shm_put
      CALL gc_ssync(nproc,info)
      CALL kk_ralltoalle(cellflux,send_map,n_items_send,sarr_len,       &
     &                    cellflux,recv_map,n_items_recv,rarr_len,      &
     &                    gid,flag,info)
      CALL gc_ssync(nproc,info)

      ndim=1 ! 1st dimension of RCELLNO
      CALL map1(send_map,recv_map,from,to,newcell,ndim,                 &
     &  n_items_send,n_items_recv)
      sarr_len=nclprc*ndim
      rarr_len=nclprc*ndim

! send out all the swapped elements of RCELLNO
      info=gc_shm_put
      CALL gc_ssync(nproc,info)
      CALL kk_ralltoalle(rcellno,send_map,n_items_send,sarr_len,        &
     &                    rcellno,recv_map,n_items_recv,rarr_len,       &
     &                    gid,flag,info)
      CALL gc_ssync(nproc,info)
      cellno=int(rcellno)

!   Fill in holes using cells at end.
      nfill=nfill+n_items_recv
      DO j=n_items_send,1,-1
        jj=hole(j)
        xx(:,jj)=xx(:,nfill)
        pos(:,jj)=pos(:,nfill)
        cellflux(:,jj)=cellflux(:,nfill)
        cellno(jj)=cellno(nfill)
        nfill=nfill-1
      END DO

      END SUBROUTINE SWAP
#endif
!   ********************************************************************
#ifdef OLDSWAP2
      SUBROUTINE SWAP2
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Swaps cell information between processors.
!-                         Internal procedure to STOCHEM
!-
!-   Inputs  :  also CELLSWAP, POS, XX, DJ, EM ,DD and DW
!-   Outputs :  also POS, XX, DJ, EM ,DD and DW
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.4    25/09/97  Created.  W.J. Collins
!  4.5    01/03/98  Reduced size of cell arrays. W.J. Collins
!  4.5    10/03/98  No longer swap Eulerian arrays. W.J. Collins
!  5.5    29/01/04  Cellflux uses num3dflux_dim dimension. M.G.Sanderson
!  5.5    11/03/04  Now calls new routine kk_ralltoalle instead of
!                   gcg_ralltoalle as former is faster. K. Ketelsen
!  6.1    21/10/04  No change.
!-
!VVV  V4.5  SWAP2 7/XII/00  LIQ included
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, DIMENSION(nclprc) :: hole
      INTEGER            :: sarr_len
      INTEGER            :: rarr_len
      INTEGER, PARAMETER :: unused = 0
      INTEGER, PARAMETER :: gid = unused
      INTEGER, PARAMETER :: flag =unused
      INTEGER :: info
      INTEGER :: n_items_send
      INTEGER :: n_items_recv
      INTEGER :: ndim
      INTEGER, DIMENSION(s_stride_in_recv_array, nclprc) :: send_map
      INTEGER, DIMENSION(r_stride_in_send_array, nclprc) :: recv_map
      INTEGER, DIMENSION(nclprc) :: from,to
      INTEGER, DIMENSION(nclprc*nproc) :: newcell
      INTEGER :: j
      INTEGER :: jj
      INTEGER :: k
      INTEGER :: swapsize
      INTEGER, DIMENSION(mxcoll) :: clswp1,nwswp1
      REAL, DIMENSION(nclprc) :: rcellno
      REAL, DIMENSION(2,nclprc*nproc) :: cellswap2
!      COMMON /CSWAP2/ CLSWP1,NWSWP1,RCELLNO

      N_ITEMS_SEND=0
      N_ITEMS_RECV=0
      NEWCELL=0
      CELLSWAP2=0
      SWAPSIZE=NCLPRC*NPROC
      RCELLNO=REAL(CELLNO)

! Sum CELLSWAP over all PEs, equivalent to broadcasting all
! non zero elements.
! Done in two sections to avoid overflow of workspace.
      DO J=1,2
        INFO=GC_SHM_PUT
        DO K=0,SWAPSIZE/MXCOLL
          CLSWP1=0
          CLSWP1(1:MIN(SWAPSIZE-K*MXCOLL,MXCOLL))=                      &
     &      CELLSWAP(J,K*MXCOLL+1:MIN((K+1)*MXCOLL,SWAPSIZE))
          CALL GC_SSYNC(NPROC,INFO)
          CALL GC_ISUM(MXCOLL,NPROC,INFO,CLSWP1)
          CALL GC_SSYNC(NPROC,INFO)
          CELLSWAP(J,K*MXCOLL+1:MIN((K+1)*MXCOLL,SWAPSIZE))=            &
     &      CLSWP1(1:MIN(SWAPSIZE-K*MXCOLL,MXCOLL))
        END DO
      END DO
! Have to remember to subtract 1 from CELLSWAP each time to get
! correct PE.
      DO J=1,NCLPRC*NPROC
        IF(CELLSWAP(1,J)-1==mype) THEN
          N_ITEMS_SEND=N_ITEMS_SEND+1
          SEND_MAP(S_DESTINATION_PE,N_ITEMS_SEND)=CELLSWAP(2,J)-1
          FROM(N_ITEMS_SEND)=J ! FROM is list of cells sent from this PE
          HOLE(N_ITEMS_SEND)=J-CELLBASE ! cell no longer on this PE
        END IF

        IF(CELLSWAP(2,J)-1==mype) THEN
          N_ITEMS_RECV=N_ITEMS_RECV+1
          RECV_MAP(R_SOURCE_PE,N_ITEMS_RECV)=CELLSWAP(1,J)-1
          TO(N_ITEMS_RECV)=J ! TO is list of cells sent to this PE
          NEWCELL(J)=NFILL+N_ITEMS_RECV ! New no. for cell on this PE
          CELLSWAP2(:,NEWCELL(J)+CELLBASE)=CELLSWAP(:,J)
        END IF
      ENDDO

! Sum NEWCELL over all PEs, equivalent to broadcasting all
! non zero elements.
      INFO=GC_SHM_PUT
      DO K=0,SWAPSIZE/MXCOLL
        NWSWP1=0
        NWSWP1(1:MIN(SWAPSIZE-K*MXCOLL,MXCOLL))=                        &
     &    NEWCELL(K*MXCOLL+1:MIN((K+1)*MXCOLL,SWAPSIZE))
        CALL GC_SSYNC(NPROC,INFO)
        CALL GC_ISUM(MXCOLL,NPROC,INFO,NWSWP1)
        CALL GC_SSYNC(NPROC,INFO)
        NEWCELL(K*MXCOLL+1:MIN((K+1)*MXCOLL,SWAPSIZE))=                 &
     &    NWSWP1(1:MIN(SWAPSIZE-K*MXCOLL,MXCOLL))
      END DO

! Get send and recv maps. Send using RALLTOALLE

      NDIM=4 ! 1st dimension of POS
      CALL MAP1(SEND_MAP,RECV_MAP,FROM,TO,NEWCELL,NDIM,                 &
     &  N_ITEMS_SEND,N_ITEMS_RECV)
      SARR_LEN=NCLPRC*NDIM
      RARR_LEN=NCLPRC*NDIM

!  send out all the swapped elements of POS
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)
      CALL KK_RALLTOALLE(POS,SEND_MAP,N_ITEMS_SEND,SARR_LEN,            &
     &                    POS,RECV_MAP,N_ITEMS_RECV,RARR_LEN,           &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)

      NDIM=NC ! 1st dimension of XX
      CALL MAP1(SEND_MAP,RECV_MAP,FROM,TO,NEWCELL,NDIM,                 &
     &  N_ITEMS_SEND,N_ITEMS_RECV)
      SARR_LEN=NCLPRC*NDIM
      RARR_LEN=NCLPRC*NDIM

! send out all the swapped elements of XX
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)

      CALL KK_RALLTOALLE(XX,SEND_MAP,N_ITEMS_SEND,SARR_LEN,             &
     &                    XX,RECV_MAP,N_ITEMS_RECV,RARR_LEN,            &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)

! send out all the swapped elements of DD
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)

      CALL KK_RALLTOALLE(DD,SEND_MAP,N_ITEMS_SEND,SARR_LEN,             &
     &                    DD,RECV_MAP,N_ITEMS_RECV,RARR_LEN,            &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)

! send out all the swapped elements of DW
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)

      CALL KK_RALLTOALLE(DW,SEND_MAP,N_ITEMS_SEND,SARR_LEN,             &
     &                    DW,RECV_MAP,N_ITEMS_RECV,RARR_LEN,            &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)

! send out all the swapped elements of EM
      NDIM=NC+2      ! 1st dimension of EM
      CALL MAP1(SEND_MAP,RECV_MAP,FROM,TO,NEWCELL,NDIM,                 &
     &  N_ITEMS_SEND,N_ITEMS_RECV)
      SARR_LEN=NCLPRC*NDIM
      RARR_LEN=NCLPRC*NDIM

      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)

      CALL KK_RALLTOALLE(EM,SEND_MAP,N_ITEMS_SEND,SARR_LEN,             &
     &                    EM,RECV_MAP,N_ITEMS_RECV,RARR_LEN,            &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)


      ndim=num3dflux_dim ! 1st dimension of CELLFLUX
      CALL MAP1(SEND_MAP,RECV_MAP,FROM,TO,NEWCELL,NDIM,                 &
     &  N_ITEMS_SEND,N_ITEMS_RECV)
      SARR_LEN=NCLPRC*NDIM
      RARR_LEN=NCLPRC*NDIM

! send out all the swapped elements of CELLFLUX
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)
      CALL KK_RALLTOALLE(CELLFLUX,SEND_MAP,N_ITEMS_SEND,SARR_LEN,       &
     &                    CELLFLUX,RECV_MAP,N_ITEMS_RECV,RARR_LEN,      &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)

      NDIM=NDJ   ! 1st dimension of DJ
      CALL MAP1(SEND_MAP,RECV_MAP,FROM,TO,NEWCELL,NDIM,                 &
     &  N_ITEMS_SEND,N_ITEMS_RECV)
      SARR_LEN=NCLPRC*NDIM
      RARR_LEN=NCLPRC*NDIM

! send out all the swapped elements of DJ
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)

      CALL KK_RALLTOALLE(DJ,SEND_MAP,N_ITEMS_SEND,SARR_LEN,             &
     &                    DJ,RECV_MAP,N_ITEMS_RECV,RARR_LEN,            &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)

      NDIM=1 ! 1st dimension of TL
      CALL MAP1(SEND_MAP,RECV_MAP,FROM,TO,NEWCELL,NDIM,                 &
     &  N_ITEMS_SEND,N_ITEMS_RECV)
      SARR_LEN=NCLPRC*NDIM
      RARR_LEN=NCLPRC*NDIM

! send out all the swapped elements of TL
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)
      CALL KK_RALLTOALLE(TL,SEND_MAP,N_ITEMS_SEND,SARR_LEN,             &
     &                    TL,RECV_MAP,N_ITEMS_RECV,RARR_LEN,            &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)

! send out all the swapped elements of QL
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)
      CALL KK_RALLTOALLE(QL,SEND_MAP,N_ITEMS_SEND,SARR_LEN,             &
     &                    QL,RECV_MAP,N_ITEMS_RECV,RARR_LEN,            &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)

! send out all the swapped elements of LIQ
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)
      CALL KK_RALLTOALLE(LIQ,SEND_MAP,N_ITEMS_SEND,SARR_LEN,            &
     &                    LIQ,RECV_MAP,N_ITEMS_RECV,RARR_LEN,           &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)

! send out all the swapped elements of M
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)
      CALL KK_RALLTOALLE(M,SEND_MAP,N_ITEMS_SEND,SARR_LEN,              &
     &                    M,RECV_MAP,N_ITEMS_RECV,RARR_LEN,             &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)

! send out all the swapped elements of RCELLNO
      INFO=GC_SHM_PUT
      CALL GC_SSYNC(NPROC,INFO)
      CALL KK_RALLTOALLE(RCELLNO,SEND_MAP,N_ITEMS_SEND,SARR_LEN,        &
     &                    RCELLNO,RECV_MAP,N_ITEMS_RECV,RARR_LEN,       &
     &                    GID,FLAG,INFO)
      CALL GC_SSYNC(NPROC,INFO)
      CELLNO=INT(RCELLNO)

!   Fill in holes using cells at end.
      NFILL=NFILL+N_ITEMS_RECV
      DO J=N_ITEMS_SEND,1,-1
        JJ=HOLE(J) ! Location of hole
        XX(:,JJ)=XX(:,NFILL)
        POS(:,JJ)=POS(:,NFILL)
        CELLFLUX(:,JJ)=CELLFLUX(:,NFILL)
        TL(JJ)=TL(NFILL)
        QL(JJ)=QL(NFILL)
        LIQ(JJ)=LIQ(NFILL)
        M(JJ)=M(NFILL)
        EM(:,JJ)=EM(:,NFILL)
        DD(:,JJ)=DD(:,NFILL)
        DW(:,JJ)=DW(:,NFILL)
        DJ(:,JJ)=DJ(:,NFILL)
        CELLNO(JJ)=CELLNO(NFILL)
        NFILL=NFILL-1
      END DO
      CELLSWAP=CELLSWAP2
      WHERE((CELLSWAP(2,:)-1) /= mype)
        CELLSWAP(1,:)=0
        CELLSWAP(2,:)=0
      END WHERE

      END SUBROUTINE SWAP2
#endif
!   ********************************************************************

!#######################################################################
      SUBROUTINE SWAP_ARRAY(data,ncols,forward)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Swap an array for load balancing
!-                         Internal procedure to STOCHEM
!-
!-   Inputs  :     data,ncols,forward
!-   Outputs :     data
!-   Controls:
!
! Current Code Owner: M.G. Sanderson
!----------------------------------------------------------------------
      USE MPL, ONLY :      &
          MPL_STATUS_SIZE, &
          MPL_REAL

      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER :: ncols
      REAL    :: data(ncols,*)
      LOGICAL :: forward

      INTEGER :: n, nreqs, nrows, lastrow, ireq(nmap), ier
      INTEGER :: istat(MPL_STATUS_SIZE,nmap)
      INTEGER :: MY_COMM

! Get communicator for MPL comms
      CALL gc_get_communicator(MY_COMM, ier)

      IF (nc_actual(mype) > nc_target(mype)) THEN

! We are sending when forward, receiving when backward

! We use non-blocking communication for this since
! we want to avoid deadlock situations due to the
! nonregular communication pattern

        lastrow = nc_target(mype) ! last row to keep
        nreqs = 0 ! number of requests
        DO n=1,nmap
          IF(balance_map(2,n) == mype) THEN
            nreqs = nreqs+1
            nrows = balance_map(3,n)
            IF (forward) THEN
              CALL MPL_Isend(data(1,lastrow+1),ncols*nrows,MPL_REAL,    &
     &                       balance_map(1,n),n,MY_COMM,                &
     &                       ireq(nreqs),ier)
            ELSE
              CALL MPL_Irecv(data(1,lastrow+1),ncols*nrows,MPL_REAL,    &
     &                       balance_map(1,n),n,MY_COMM,                &
     &                       ireq(nreqs),ier)
            END IF
            lastrow = lastrow+nrows
          END IF
        END DO

! Wait for all communication to complete

        CALL MPL_Waitall(nreqs,ireq,istat,ier)

      ELSE

! We are receiving when forward, sending when backward

! For avoiding deadlock situations it is sufficient when
! one partner uses non-blocking communication, so there
! is no need to use it here, too.

        lastrow = nc_actual(mype) ! last row on this PE
        DO n=1,nmap
          IF(balance_map(1,n) == mype) THEN
            nrows = balance_map(3,n)
            IF (forward) THEN
              CALL MPL_Recv(data(1,lastrow+1),ncols*nrows,MPL_REAL,     &
     &                      balance_map(2,n),n,MY_COMM,istat,ier)
            ELSE
              CALL MPL_Send(data(1,lastrow+1),ncols*nrows,MPL_REAL,     &
     &                      balance_map(2,n),n,MY_COMM,ier)
            END IF
            lastrow = lastrow+nrows
          END IF
        END DO
      END IF

      END SUBROUTINE SWAP_ARRAY
      SUBROUTINE ZERO_STATS
!-----------------------------------------------------------------------
!
!     Zero output statistics.
!     Internal procedure to STOCHEM.
!
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.1    14/03/00  Created.  D.S. Stevenson
!  6.1    21/10/04  No change.
!
!VVV  V2.3  zero_stats 14/III/00
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      tflux=0.
      totflu=0.
      totavg=0.
      period=0.
      navg=0
      mconc=0.
      sdconc=0.
      np=0
      nph=0
      nm=0
      nstat=0
! DEPENDS ON: mcalc
      CALL MCALC(M0,totm0,xx,clist,ipos,pos,tropz,nchem,nfill)

      END SUBROUTINE ZERO_STATS
      END SUBROUTINE STOCHEM
#endif
