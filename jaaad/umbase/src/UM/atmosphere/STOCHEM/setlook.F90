#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Subroutine to set up the pp lookup table.
!
!     Subroutine Interface:
      SUBROUTINE SETLOOK(month,year,daym0,daymas,flist,clist,stash_code,&
     &  nchem,ndim,flheader,pplook)
!
!
!      Implicit None
!
! Description:
!   This routine takes the ozone field supplied in the ancillary
!
! Method:
!   This routine sets up the lookup table depending on the
!   contents of the CLIST and FLIST arrays.
!
! Current Code Owner: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.5    03/08/99  Created.  C.E. Johnson
!  5.5    27/01/04  Changed STASH sections to avoid clash with river
!                   routing diagnostics. M.G. Sanderson
!  6.1    08/03/04  Original Code adapted from earlier STOCHEM vns.
!                   C.E. Johnson
!  6.2    01/03/06  Now outputs STASH code from strings in module.
!                   Also sets lblev values.   M.G. Sanderson.
!
! Code description:
!   FORTRAN 90
!   This code is written to the programming standards of version 6
!   of UMDP3.
!
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_OUT
      IMPLICIT NONE
!-----------------------------------------------------------------------
#include "typsize.h"
      INTEGER,                       INTENT(IN)         :: month
      INTEGER,                       INTENT(IN)         :: year
      INTEGER,                       INTENT(IN)         :: ndim
      INTEGER,                       INTENT(IN)         :: nchem
      INTEGER, DIMENSION(2,numflux), INTENT(IN)         :: flist
      INTEGER, DIMENSION(numchem),   INTENT(IN)         :: clist
      INTEGER, DIMENSION(numflux),   INTENT(IN)         :: stash_code
      REAL,                          INTENT(IN)         :: daym0
      REAL,                          INTENT(IN)         :: daymas

      INTEGER, DIMENSION(len_fixhd), INTENT(IN)         :: flheader
      INTEGER, DIMENSION(len1_lookup,ndim), INTENT(OUT) :: pplook

      INTEGER  :: icode
      CHARACTER*80 :: cmessage

      INTEGER                  :: i
      INTEGER                  :: j
      INTEGER                  :: l
      INTEGER                  :: ll
      INTEGER, DIMENSION(ndim) :: ilev
      INTEGER :: Stash_Section_1=25  ! Concentrations
      INTEGER :: Stash_Section_2=24  ! Cell positions and 3D fluxes (was
!  Changed to 24 as section 26 is river routing diagnostics in UMUI
      INTEGER                  :: item
      INTEGER                  :: exptcode
      INTEGER, DIMENSION(ndim) :: fc
      INTEGER, DIMENSION(ndim) :: stash

      REAL                     :: deltax
      REAL                     :: deltay
      CHARACTER(LEN=14)        :: filename
      CHARACTER(LEN=5)         :: run_id

      INTEGER, DIMENSION((numchem+max3dflux+1)*nlev) :: ippval
      REAL, DIMENSION((numchem+max3dflux+1)*nlev)    :: rppval
      EQUIVALENCE (ippval,rppval)

#include "cntl_io.h"
#include "c_mdi.h"
#include "clookadd.h"
      EXTERNAL EXPT_ENC

! Set field codes and STASH codes by interrogating FLIST. For fluxes,
! field codes are consecutive and do not depend on the flux.
! Cell position array is given field code 1800 (for consistency
! with older versions of STOCHEM) and STASH code 24001. Then, all
! 3D fluxes have the item/field code starting at 2, i.e. first
! 3D flux has field code = 1802, STASH code 24002, and so on.
!
      item = 1                    ! Item must be >= 1
      fc(1:nlev) = 1800           ! Cell position array
      stash(1:nlev) = 24401       ! Hard-code this for now

      j = 1                       ! Counter for fields output to STASH
      DO l=1,numflux              ! 3-D fluxes
        IF (flist(2,l) > 0) THEN
          j = j + 1
          item = j
          fc(j*nlev-(nlev-1):j*nlev) = 1800 + item
          stash(j*nlev-(nlev-1):j*nlev) = stash_code(l)
        END IF
      END DO

      DO ll=1,nchem               ! Concentration fields
        l = clist(ll)
        j = j + 1
        fc(j*nlev-(nlev-1):j*nlev) = fcodes(l)
        item = fcodes(l) - 1700
        stash(j*nlev-(nlev-1):j*nlev) = stash_section_1*1000 + item
      END DO

      IF (j*nlev > ndim) THEN
        cmessage = 'Too many output fluxes requested in flist'
        WRITE(6,*) cmessage,j,nlev,ndim
! DEPENDS ON: ereport
        CALL EREPORT('SETLOOK',1,cmessage)
      END IF

      icode = 0
      pplook = 0
      run_id(1:4) = expt_id
      run_id(5:5) = job_id
! DEPENDS ON: expt_enc
      CALL EXPT_ENC(run_id,exptcode,icode,cmessage)
      IF (icode /= 0) THEN
! DEPENDS ON: ereport
        CALL EREPORT('EXPT_ENC',1,cmessage)
      END IF

      pplook(lbyr,:)=   year                      ! | First
      pplook(lbmon,:)=  month                     ! | Validity
      pplook(lbdat,:)=  INT(daym0)                ! | Time
      pplook(lbhr,:)=   INT((daym0-INT(daym0))*24)
      pplook(lbmin,:)=  MOD(INT(daym0*24*60),60)
      pplook(lbyrd,:)=  year                      ! + Last
      pplook(lbmond,:)= month                     ! + Validity
      pplook(lbdatd,:)= INT(daymas)               ! + Time
      pplook(lbhrd,:)=  INT((daymas-INT(daymas))*24)
      pplook(lbmind,:)= MOD(INT(daymas*24*60),60)
      pplook(lbtim,:)=  3*100+2*10+2              ! time indicator
      pplook(lblrec,:)= mnlat*nlong               ! Data words
      pplook(lbcode,:)= 1                         ! Regular Lat x long g
      pplook(lbhem,:)=  0                         ! Global field
      pplook(lbrow,:)=  mnlat                     ! No. of rows
      pplook(lbnpt,:)=  nlong                     ! No. of cols
      pplook(lbpack,:)= 02000                     ! No packing
      pplook(lbrel,:)=  2                         ! Header release numbe
      pplook(lbfc,:)=   fc                        ! Field code
      pplook(lbproc,:)= 128                       ! Time mean
      pplook(lbvc,:)=   65                        ! Hybrid height grid
      pplook(lbexp,:)=  exptcode                  ! Expt. id.
      pplook(lbnrec,:)= ((pplook(lblrec,:)+um_sector_size-1)/           &
     &  um_sector_size)*um_sector_size            ! Rounded length
      pplook(lbegin,1)=flheader(150)-1+64*ndim
      pplook(lbegin,1)=((pplook(lbegin,1)+um_sector_size-1)/            &
     &  um_sector_size)*um_sector_size            ! Rounded length

      DO i=2,ndim
        pplook(lbegin,i)=pplook(lbegin,i-1)+pplook(lbnrec,i)
      END DO

      DO i = 1, ndim
        pplook(lblev,i) = 1 + MOD(i-1, nlev)      ! Vert level number
      END DO
      pplook(data_type,:)=1                  ! Real data
      pplook(naddr,:)=    pplook(lbegin,:)   ! Start address of data
      pplook(item_code,:)=stash              ! Stash code

      ilev = MOD((/(i-1,i=1,ndim)/),nlev)+1
      rppval(1:ndim)=Zsea_stochem(ilev)      ! Zsea value of upper bound
      pplook(bulev,:)=    ippval(1:ndim)
      rppval(1:ndim)=Ck_stochem(ilev)        ! Ck value of upper boundar
      pplook(bhulev,:)=   ippval(1:ndim)
      rppval(1:ndim)=Zsea_stochem_half(ilev) ! Zsea value of level
      pplook(blev,:)=     ippval(1:ndim)
      rppval(1:ndim)=Zsea_stochem(ilev-1)    ! Zsea value of lower bound
      pplook(brlev,:)=    ippval(1:ndim)
      rppval(1:ndim)=Ck_stochem_half(ilev)   ! Ck value of level
      pplook(bhlev,:)=    ippval(1:ndim)
      rppval(1:ndim)=Ck_stochem(ilev-1)      ! Ck value of lower boundar
      pplook(bhrlev,:)=   ippval(1:ndim)
      rppval(1:ndim)=90.                     ! BPLAT
      pplook(bplat,:)=    ippval(1:ndim)
      rppval(1:ndim)=0.
      pplook(bplon,:)=    ippval(1:ndim)     ! BPLON
      deltay=180.0/REAL(mnlat)               ! latitudinal interval
      rppval(1:ndim)=-90.0-deltay/2.0        ! zeroth latitude
      pplook(bzy,:)=      ippval(1:ndim)
      rppval(1:ndim)=deltay                  ! Latitudinal interval
      pplook(bdy,:)=      ippval(1:ndim)
      deltax=360.0/REAL(nlong)               ! longtitudinal interval
      rppval(1:ndim)=0.0-deltax/2.0          ! Zeroth longtitude
      pplook(bzx,:)=      ippval(1:ndim)
      rppval(1:ndim)=deltax                  ! Longtitudinal interval
      pplook(bdx,:)=      ippval(1:ndim)
      rppval(1:ndim)=rmdi                    ! Missing data indicator
      pplook(bmdi,:)=     ippval(1:ndim)
      rppval(1:ndim)=1.0                     ! MKS scaling factor
      pplook(bmks,:)=     ippval(1:ndim)

      END SUBROUTINE SETLOOK
#endif
