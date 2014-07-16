#if defined(A34_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module defining ASAD arrays, variables, and parameters
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Current Code Owner:       Colin Johnson
!
!  Code Description:
!   Language:  FORTRAN 90 
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------

MODULE ASAD_MOD
IMPLICIT NONE
SAVE

! Interface section
INTERFACE asad_mod_init
  MODULE PROCEDURE asad_mod_init
END INTERFACE asad_mod_init

INTERFACE asad_mod_final
  MODULE PROCEDURE asad_mod_final
END INTERFACE asad_mod_final

REAL, ALLOCATABLE :: wp(:)           ! water vapour field
REAL, ALLOCATABLE :: dpd(:,:)
REAL, ALLOCATABLE :: dpw(:,:)
REAL, ALLOCATABLE :: emr(:,:)
REAL, ALLOCATABLE :: fj(:,:,:)       ! Full jacobian
REAL, ALLOCATABLE :: qa(:,:)
REAL, ALLOCATABLE :: ratio(:,:)
REAL, ALLOCATABLE :: p(:)
REAL, ALLOCATABLE :: t(:)
REAL, ALLOCATABLE :: t300(:)
REAL, ALLOCATABLE :: tnd(:)
REAL, ALLOCATABLE :: pmintnd(:)
REAL, ALLOCATABLE :: f(:,:)
REAL, ALLOCATABLE :: fdot(:,:)
REAL, TARGET, ALLOCATABLE :: pd(:,:)
REAL, POINTER     :: prod(:,:)
REAL, POINTER     :: slos(:,:)
REAL, ALLOCATABLE :: y(:,:)
REAL, ALLOCATABLE :: ydot(:,:)
REAL, ALLOCATABLE :: ftilde(:,:)     ! lower order solution 
REAL, ALLOCATABLE :: f_no(:,:)       ! partioning factors
REAL, ALLOCATABLE :: f_no2(:,:)      !    between
REAL, ALLOCATABLE :: f_no3(:,:)      ! NO, NO2, and NO3
REAL, ALLOCATABLE :: ej(:,:)
REAL, ALLOCATABLE :: rk(:,:)
REAL, ALLOCATABLE :: prk(:,:)
REAL, ALLOCATABLE :: deriv(:,:,:)
REAL, ALLOCATABLE :: sc(:,:,:)       ! Species concentration array
REAL, ALLOCATABLE :: za(:)           ! Solar zenith angle
REAL, ALLOCATABLE :: co3(:)          ! Column ozone
REAL, ALLOCATABLE :: lati(:)         ! Latitude
REAL, ALLOCATABLE :: sphno3(:)       ! Amount of HNO3 in solid phase
REAL, ALLOCATABLE :: sph2o(:)        ! Amount of H2O in solid phase
REAL, ALLOCATABLE :: depvel(:,:,:)
REAL, ALLOCATABLE :: k298(:)
REAL, ALLOCATABLE :: dhr(:)
REAL, ALLOCATABLE :: kd298(:)
REAL, ALLOCATABLE :: ddhr(:)
REAL, ALLOCATABLE :: ab(:,:)
REAL, ALLOCATABLE :: at(:,:)
REAL, ALLOCATABLE :: aj(:,:)
REAL, ALLOCATABLE :: ah(:,:)
REAL, ALLOCATABLE :: frpb(:)         ! fractional product array
REAL, ALLOCATABLE :: ztabpd(:,:)
REAL, ALLOCATABLE :: shno3(:)        ! No. density type 1 psc solid phase hno3
REAL, ALLOCATABLE :: sh2o(:)         ! No. density type 2 psc solid phase h2o
REAL, ALLOCATABLE :: fpsc1(:)        ! 1.0 if type 1 psc's are present, else 0
REAL, ALLOCATABLE :: fpsc2(:)        ! 1.0 if type 2 psc's are present, else 0
REAL, ALLOCATABLE :: spfj(:,:)       ! Sparse full Jacobian

INTEGER, ALLOCATABLE :: madvtr(:)
INTEGER, ALLOCATABLE :: majors(:)
INTEGER, ALLOCATABLE :: moffam(:)
INTEGER, ALLOCATABLE :: nodd(:)
INTEGER, ALLOCATABLE :: nltrf(:)
INTEGER, ALLOCATABLE :: nltr3(:)
INTEGER, ALLOCATABLE :: ipa(:,:)     ! Pivot information for solving jacobian
INTEGER, ALLOCATABLE :: ipa2(:)
INTEGER, ALLOCATABLE :: nltrim(:,:)
INTEGER, ALLOCATABLE :: nlpdv(:,:)
INTEGER, ALLOCATABLE :: nfrpx(:)     ! index to product array
INTEGER, ALLOCATABLE :: ntabfp(:,:)
INTEGER, ALLOCATABLE :: ntabpd(:,:)
INTEGER, ALLOCATABLE :: npdfr(:,:)
INTEGER, ALLOCATABLE :: ngrp(:,:)
INTEGER, ALLOCATABLE :: njcgrp(:,:)
INTEGER, ALLOCATABLE :: nprdx3(:,:,:)
INTEGER, ALLOCATABLE :: nprdx2(:,:)
INTEGER, ALLOCATABLE :: nprdx1(:)
INTEGER, ALLOCATABLE :: njacx3(:,:,:)
INTEGER, ALLOCATABLE :: njacx2(:,:)
INTEGER, ALLOCATABLE :: njacx1(:)
INTEGER, ALLOCATABLE :: nmpjac(:)
INTEGER, ALLOCATABLE :: npjac1(:,:)
INTEGER, ALLOCATABLE :: nbrkx(:)
INTEGER, ALLOCATABLE :: ntrkx(:)
INTEGER, ALLOCATABLE :: nprkx(:)
INTEGER, ALLOCATABLE :: nhrkx(:)
INTEGER, ALLOCATABLE :: nlall(:)
INTEGER, ALLOCATABLE :: nlstst(:)
INTEGER, ALLOCATABLE :: nlf(:)
INTEGER, ALLOCATABLE :: nlmajmin(:)
INTEGER, ALLOCATABLE :: nldepd(:)
INTEGER, ALLOCATABLE :: nldepw(:)
INTEGER, ALLOCATABLE :: nlemit(:)
INTEGER, ALLOCATABLE :: nldepx(:)
INTEGER, ALLOCATABLE :: njcoth(:,:)
INTEGER, ALLOCATABLE :: nmzjac(:)
INTEGER, ALLOCATABLE :: nzjac1(:,:)
INTEGER, ALLOCATABLE :: njcoss(:,:)
INTEGER, ALLOCATABLE :: nmsjac(:)
INTEGER, ALLOCATABLE :: nsjac1(:,:)
INTEGER, ALLOCATABLE :: nsspt(:)
INTEGER, ALLOCATABLE :: nspi(:,:)
INTEGER, ALLOCATABLE :: nsspi(:,:)
INTEGER, ALLOCATABLE :: nssi(:)
INTEGER, ALLOCATABLE :: nssrt(:)
INTEGER, ALLOCATABLE :: nssri(:,:)
INTEGER, ALLOCATABLE :: nssrx(:,:)
! sparse algebra
INTEGER, ALLOCATABLE :: pointer(:,:)  ! Map of nonzero entries
INTEGER, ALLOCATABLE :: pointer1(:,:) ! modified map (after decomposition)
INTEGER, ALLOCATABLE :: pointer2(:,:) ! Map of nonzero entries, before reordering
INTEGER, ALLOCATABLE :: ro(:)         ! reordering of tracers to minimize fill-in

REAL, ALLOCATABLE :: zonal_total(:,:,:)   ! allocated in ukca_budget2

INTEGER, ALLOCATABLE :: ilcf(:)
INTEGER, ALLOCATABLE :: ilss(:)
INTEGER, ALLOCATABLE :: ilct(:)
INTEGER, ALLOCATABLE :: ilftr(:)
INTEGER, ALLOCATABLE :: ilft(:)
INTEGER, ALLOCATABLE :: ilstmin(:)
! budget2 vars:
INTEGER, ALLOCATABLE :: rpartner_bimol(:,:)
INTEGER, ALLOCATABLE :: rpartner_termol(:,:)
INTEGER, ALLOCATABLE :: rpartner_hetero(:,:)
INTEGER, ALLOCATABLE :: rpartner_wd(:)
INTEGER, ALLOCATABLE :: rpartner_dd(:)
INTEGER, ALLOCATABLE :: rpartner_phot(:)

LOGICAL, ALLOCATABLE :: linfam(:,:)
LOGICAL, ALLOCATABLE :: ldepd(:)     ! T for dry deposition
LOGICAL, ALLOCATABLE :: ldepw(:)     ! T for wet deposition
LOGICAL, ALLOCATABLE :: lemit(:)     ! T for emission

CHARACTER(len=10), ALLOCATABLE :: advt(:)
CHARACTER(len=10), ALLOCATABLE :: family(:)
CHARACTER(len=10), ALLOCATABLE :: speci(:)
CHARACTER(len=2),  ALLOCATABLE :: ctype(:)
CHARACTER(len=10), ALLOCATABLE :: spb(:,:)
CHARACTER(len=10), ALLOCATABLE :: spt(:,:)
CHARACTER(len=10), ALLOCATABLE :: spj(:,:)
CHARACTER(len=10), ALLOCATABLE :: sph(:,:)

REAL, PARAMETER    :: pmin = 1.0e-20
REAL, PARAMETER    :: ptol = 1.0e-5
REAL, PARAMETER    :: ftol = 1.0e-3

INTEGER, PARAMETER :: kfphot=0
INTEGER, PARAMETER :: jpss = 13
INTEGER, PARAMETER :: jpssr = 31
INTEGER, PARAMETER :: nvar = 17
INTEGER, PARAMETER :: nllv = 17
INTEGER, PARAMETER :: ninv = 23
INTEGER, PARAMETER :: nout=71
!     nout:        Fortran channel for output in subroutine OUTVMR
INTEGER, PARAMETER :: jpem = 9      ! IS THIS A CONSTANT/USED ??
INTEGER, PARAMETER :: jddept = 6
!     jddept:      Number of time periods used in dry deposition i.e.
!                  summer(day,night,24h ave), winter(day,night,24h ave)
INTEGER, PARAMETER :: jddepc = 5
!     jddepc:      Number of land use categories used in dry dep.
INTEGER, PARAMETER :: jpdwio = 56
!     jpdwio       Fortran i/o unit to read/write anything to do with
!                  wet/dry deposition
INTEGER, PARAMETER :: jpemio = 57
!     jpemio       Fortran i/o unit to read in emissions
INTEGER, PARAMETER :: jpfrpd=100
INTEGER, PARAMETER :: jpab    = 3
INTEGER, PARAMETER :: jpat    = 7
INTEGER, PARAMETER :: jpaj    = 3
INTEGER, PARAMETER :: jpah    = 3
INTEGER, PARAMETER :: jpspb   = 6
INTEGER, PARAMETER :: jpspt   = 4
INTEGER, PARAMETER :: jpspj   = 6
INTEGER, PARAMETER :: jpsph   = 6
INTEGER, PARAMETER :: jpmsp   = jpspb
INTEGER, PARAMETER :: jppjac  = 10
INTEGER, PARAMETER :: jpkargs = 10
INTEGER, PARAMETER :: jprargs = 10
INTEGER, PARAMETER :: jpcargs = 1
INTEGER, PARAMETER :: jpfrpb  = (jpspb-2)*30
INTEGER, PARAMETER :: jpcio   = 55
INTEGER, PARAMETER :: spfjsize_max = 700 ! maximum number of
                                         ! nonzero matrix elements


CHARACTER(LEN=2),  PARAMETER :: jpfm = 'FM'
CHARACTER(LEN=2),  PARAMETER :: jpif = 'FT'
CHARACTER(LEN=2),  PARAMETER :: jpsp = 'TR'
CHARACTER(LEN=2),  PARAMETER :: jpna = 'SS'
CHARACTER(LEN=2),  PARAMETER :: jpco = 'CT'
CHARACTER(LEN=2),  PARAMETER :: jpcf = 'CF'
CHARACTER(LEN=15), PARAMETER :: FDEPVEL = 'depvel.d'
CHARACTER(LEN=15), PARAMETER :: FHENRY  = 'henry.d'
CHARACTER(LEN=15), PARAMETER :: CHCH = 'chch.d'
CHARACTER(LEN=15), PARAMETER :: RATB = 'ratb.d'
CHARACTER(LEN=15), PARAMETER :: RATT = 'ratt.d'
CHARACTER(LEN=15), PARAMETER :: RATJ = 'ratj.d'
CHARACTER(LEN=15), PARAMETER :: RATH = 'rath.d'

LOGICAL, PARAMETER :: lvmr=.true.    ! T for volume mixing ratio

REAL :: fch4,fco2,fh2,fn2,fo2
REAL    :: cdt                       ! chemistry timestep
REAL    :: peps                      ! 
REAL    :: dtime                     !

INTEGER :: kcdt
INTEGER :: nrsteps
INTEGER :: nitnr
INTEGER :: nitfg
INTEGER :: ntrf, ntr3
INTEGER :: jlp            ! level to be integrated
INTEGER :: nuni, nsst
INTEGER :: ncsteps
INTEGER :: nit0=20        ! ftoy iterations with method=0
INTEGER :: nfphot
INTEGER :: jsubs
INTEGER :: jlst
INTEGER :: method         ! chemistry integration method
INTEGER :: nnfrp
INTEGER :: nfrpd
INTEGER :: nstst, nf, ndepd, ndepw, nemit
INTEGER :: ntro3, ntroh, ntrho2, ntrno 
INTEGER :: nspo1d, nspo3p, nspo3, nspoh
INTEGER :: nspho2, nspno, nspn, nsph

LOGICAL :: lsvjac         ! Flag for saving jacobian if not recalculated
LOGICAL :: ljacx

CONTAINS

! ######################################################################
SUBROUTINE asad_mod_init

! To allocate and initialise ASAD arrays and variables

IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"

ALLOCATE(wp(theta_field_size))
ALLOCATE(linfam(theta_field_size,0:jpctr))
ALLOCATE(madvtr(jpspec))
ALLOCATE(majors(jpctr))
ALLOCATE(moffam(jpspec))
ALLOCATE(nodd(jpspec))
ALLOCATE(nltrf(jpctr))
ALLOCATE(nltr3(jpctr))
ALLOCATE(advt(jpctr))
ALLOCATE(family(jpspec))
ALLOCATE(speci(jpspec))
ALLOCATE(ctype(jpspec))
ALLOCATE(dpd(theta_field_size,jpspec))
ALLOCATE(dpw(theta_field_size,jpspec))
ALLOCATE(emr(theta_field_size,jpspec))
ALLOCATE(ipa(theta_field_size,jpctr))
ALLOCATE(ipa2(jpctr))
ALLOCATE(fj(theta_field_size,jpctr,jpctr))
ALLOCATE(qa(theta_field_size,jpspec))
ALLOCATE(ratio(theta_field_size,jpspec))
ALLOCATE(p(theta_field_size))
ALLOCATE(t(theta_field_size))
ALLOCATE(t300(theta_field_size))
ALLOCATE(tnd(theta_field_size))
ALLOCATE(pmintnd(theta_field_size))
ALLOCATE(f(theta_field_size,jpctr))
ALLOCATE(fdot(theta_field_size,jpctr))
ALLOCATE(pd(theta_field_size,2*jpspec))
ALLOCATE(y(theta_field_size,jpspec))
ALLOCATE(ydot(theta_field_size,jpspec))
ALLOCATE(ftilde(theta_field_size, jpctr))
ALLOCATE(f_no(theta_field_size,model_levels))
ALLOCATE(f_no2(theta_field_size,model_levels))
ALLOCATE(f_no3(theta_field_size,model_levels))          
ALLOCATE(ej(theta_field_size,jpctr))
ALLOCATE(rk(theta_field_size,jpnr))
ALLOCATE(prk(theta_field_size,jpnr))
ALLOCATE(deriv(theta_field_size,4,4))
ALLOCATE(nspi(jpnr,jpmsp))
ALLOCATE(nsspt(jpss))
ALLOCATE(nsspi(jpss,jpssr))
ALLOCATE(nssi(jpss))
ALLOCATE(nssrt(jpss))
ALLOCATE(nssri(jpss,jpssr))
ALLOCATE(nssrx(jpss,jpssr)) 
ALLOCATE(sc(theta_field_size,model_levels,jpspec))
ALLOCATE(za(theta_field_size))
ALLOCATE(co3(theta_field_size))
ALLOCATE(lati(theta_field_size))
ALLOCATE(sphno3(theta_field_size))
ALLOCATE(sph2o(theta_field_size))
ALLOCATE(depvel(jddept,jddepc,jpdd))
ALLOCATE(k298(jpdw))
ALLOCATE(dhr(jpdw))
ALLOCATE(kd298(jpdw))
ALLOCATE(ddhr(jpdw))
ALLOCATE(ldepd(jpspec))
ALLOCATE(ldepw(jpspec))
ALLOCATE(lemit(jpspec))
ALLOCATE(nltrim(0:jpctr,3))
ALLOCATE(nlpdv((jpspj-2)*jppj,2))
ALLOCATE(ab(jpbk+1,jpab))
ALLOCATE(at(jptk+1,jpat))
ALLOCATE(aj(jppj+1,jpaj))
ALLOCATE(ah(jphk+1,jpah))
ALLOCATE(spb(jpbk+1,jpspb))
ALLOCATE(spt(jptk+1,jpspt))
ALLOCATE(spj(jppj+1,jpspj))
ALLOCATE(sph(jphk+1,jpsph))
ALLOCATE(frpb(jpfrpb))
ALLOCATE(ztabpd(jpfrpd,2))
ALLOCATE(nfrpx(jpnr))
ALLOCATE(ntabfp(jpfrpb,3))
ALLOCATE(ntabpd(jpfrpd,3))
ALLOCATE(npdfr(jpnr,2))
ALLOCATE(ngrp(2*jpspec,3))
ALLOCATE(njcgrp(jpctr,3))
ALLOCATE(nprdx3(3,(jpnr/(3*3))+3*3,2*jpspec))
ALLOCATE(nprdx2(2,2*jpspec))
ALLOCATE(nprdx1(2*jpspec))
ALLOCATE(njacx3(3,(jpnr/(3*3))+3*3,jpctr))
ALLOCATE(njacx2(2,jpctr))
ALLOCATE(njacx1(jpctr))
ALLOCATE(nmpjac(jpctr))
ALLOCATE(npjac1(jppjac,jpctr))
ALLOCATE(nbrkx(jpbk+1))
ALLOCATE(ntrkx(jptk+1))
ALLOCATE(nprkx(jppj+1))
ALLOCATE(nhrkx(jphk+1))
ALLOCATE(nlall(jpspec))
ALLOCATE(nlstst(jpspec))
ALLOCATE(nlf(jpspec))
ALLOCATE(nlmajmin(jpspec))
ALLOCATE(nldepd(jpspec))
ALLOCATE(nldepw(jpspec))
ALLOCATE(nlemit(jpspec))
ALLOCATE(nldepx(jpspec))
ALLOCATE(njcoth(jpnr,jpmsp))
ALLOCATE(nmzjac(jpctr))
ALLOCATE(nzjac1(jpnr,jpctr))
ALLOCATE(njcoss(jpnr,jpmsp))
ALLOCATE(nmsjac(jpctr))
ALLOCATE(nsjac1(jpnr,jpctr))
ALLOCATE(shno3(theta_field_size))
ALLOCATE(sh2o(theta_field_size))
ALLOCATE(fpsc1(theta_field_size))
ALLOCATE(fpsc2(theta_field_size))

! the following had save attribs.
ALLOCATE(ilcf(jpspec))
ALLOCATE(ilss(jpspec))
ALLOCATE(ilct(jpspec))
ALLOCATE(ilftr(jpspec))
ALLOCATE(ilft(jpspec))
ALLOCATE(ilstmin(jpspec))

! EQUIVALENCE ( pd(1,1), prod(1,1) )
! EQUIVALENCE ( pd(1,jpspec+1), slos(1,1) )
prod => pd(:,1:jpspec)
slos => pd(:,jpspec+1:2*jpspec)

! Set integration method
IF (L_ukca_strat .OR. L_ukca_strattrop .OR. L_ukca_stratcfc) THEN
  method = 3
ELSE
  method = 1
ENDIF

! Initialize variables that may be changed in cinit (was in ini_runcnsts)
kcdt    = 900
nrsteps = 15
nitnr   = 10
nitfg   = 10
dtime   = 900.0

IF (method == 3) THEN
  ALLOCATE(pointer(jpctr, jpctr))
  ALLOCATE(pointer1(jpctr, jpctr))
  ALLOCATE(pointer2(jpctr, jpctr))
  ALLOCATE(ro(jpctr))
  ALLOCATE(spfj(theta_field_size,spfjsize_max))
ENDIF

IF (L_ukca_budget2) THEN         ! vars for budget2
  ALLOCATE(rpartner_bimol(jpbk,2))
  ALLOCATE(rpartner_termol(jptk,2))
! set to zero to cope with jphk=0 condition
  ALLOCATE(rpartner_hetero(0:jphk,2))
  ALLOCATE(rpartner_wd(jpdw))
  ALLOCATE(rpartner_dd(jpdd))
  ALLOCATE(rpartner_phot(jppj))
ENDIF

RETURN
END SUBROUTINE asad_mod_init

! ######################################################################
SUBROUTINE asad_mod_final

! To deallocate ASAD arrays

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

DEALLOCATE(wp)
DEALLOCATE(linfam)
DEALLOCATE(madvtr)
DEALLOCATE(majors)
DEALLOCATE(moffam)
DEALLOCATE(nodd)
DEALLOCATE(nltrf)
DEALLOCATE(nltr3)
DEALLOCATE(advt)
DEALLOCATE(family)
DEALLOCATE(speci)
DEALLOCATE(ctype)
DEALLOCATE(dpd)
DEALLOCATE(dpw)
DEALLOCATE(emr)
DEALLOCATE(ipa)
DEALLOCATE(ipa2)
DEALLOCATE(fj)
DEALLOCATE(qa)
DEALLOCATE(ratio)
DEALLOCATE(p)
DEALLOCATE(t)
DEALLOCATE(t300)
DEALLOCATE(tnd)
DEALLOCATE(pmintnd)
DEALLOCATE(f)
DEALLOCATE(fdot)
DEALLOCATE(pd)
DEALLOCATE(y)
DEALLOCATE(ydot)
DEALLOCATE(ftilde)
DEALLOCATE(f_no)
DEALLOCATE(f_no2)
DEALLOCATE(f_no3)
DEALLOCATE(ej)
DEALLOCATE(rk)
DEALLOCATE(prk)
DEALLOCATE(deriv)
DEALLOCATE(nspi)
DEALLOCATE(nsspt)
DEALLOCATE(nsspi)
DEALLOCATE(nssi)
DEALLOCATE(nssrt)
DEALLOCATE(nssri)
DEALLOCATE(nssrx)
DEALLOCATE(sc)
DEALLOCATE(za)
DEALLOCATE(co3)
DEALLOCATE(lati)
DEALLOCATE(sphno3)
DEALLOCATE(sph2o)
DEALLOCATE(depvel)
DEALLOCATE(k298)
DEALLOCATE(dhr)
DEALLOCATE(kd298)
DEALLOCATE(ddhr)
DEALLOCATE(ldepd)
DEALLOCATE(ldepw)
DEALLOCATE(lemit)
DEALLOCATE(nltrim)
DEALLOCATE(nlpdv)
DEALLOCATE(ab)
DEALLOCATE(at)
DEALLOCATE(aj)
DEALLOCATE(ah)
DEALLOCATE(spb)
DEALLOCATE(spt)
DEALLOCATE(spj)
DEALLOCATE(sph)
DEALLOCATE(frpb)
DEALLOCATE(ztabpd)
DEALLOCATE(nfrpx)
DEALLOCATE(ntabfp)
DEALLOCATE(ntabpd)
DEALLOCATE(npdfr)
DEALLOCATE(ngrp)
DEALLOCATE(njcgrp)
DEALLOCATE(nprdx3)
DEALLOCATE(nprdx2)
DEALLOCATE(nprdx1)
DEALLOCATE(njacx3)
DEALLOCATE(njacx2)
DEALLOCATE(njacx1)
DEALLOCATE(nmpjac)
DEALLOCATE(npjac1)
DEALLOCATE(nbrkx)
DEALLOCATE(ntrkx)
DEALLOCATE(nprkx)
DEALLOCATE(nhrkx)
DEALLOCATE(nlall)
DEALLOCATE(nlstst)
DEALLOCATE(nlf)
DEALLOCATE(nlmajmin)
DEALLOCATE(nldepd)
DEALLOCATE(nldepw)
DEALLOCATE(nlemit)
DEALLOCATE(nldepx)
DEALLOCATE(njcoth)
DEALLOCATE(nmzjac)
DEALLOCATE(nzjac1)
DEALLOCATE(njcoss)
DEALLOCATE(nmsjac)
DEALLOCATE(nsjac1)
DEALLOCATE(shno3)
DEALLOCATE(sh2o)
DEALLOCATE(fpsc1)
DEALLOCATE(fpsc2)

IF (method == 3) THEN       ! sparse vars
  DEALLOCATE(pointer)
  DEALLOCATE(pointer1)
  DEALLOCATE(pointer2)
  DEALLOCATE(ro)
  DEALLOCATE(spfj)
ENDIF

IF (L_ukca_budget2) THEN         ! vars for budget2
  DEALLOCATE(rpartner_bimol)
  DEALLOCATE(rpartner_termol)
  DEALLOCATE(rpartner_hetero)
  DEALLOCATE(rpartner_wd)
  DEALLOCATE(rpartner_dd)
  DEALLOCATE(rpartner_phot)
ENDIF

RETURN
END SUBROUTINE asad_mod_final

END MODULE
#endif
