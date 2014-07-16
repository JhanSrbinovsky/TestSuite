#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EGET(estore,em,emiss,ipos,pos,bl,orog,nnn,nbl,month,   &
     & time,so2em,missings,be7em,missingbe7,be10em,                     &
     & missingbe10,lnoxem,missinglnox,acnoxem,missingacnox,nfill)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : GET EMISSIONS
!-
!-   Inputs  : EMISS,ESTORE,NNN,MISSING
!-   Outputs : EM
!-   Controls:
!-
!
! Current Code Owner: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.1    06/02/96  Created.  W.J. Collins
!  5.5    29/11/01  SO2EM, MISSINGS added. C.E. Johnson
!  5.5    13/02/04  Vectorised code. K. Ketelsen
!  6.1    22/10/04  No change.
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER,                            INTENT(IN) :: month
      INTEGER,                            INTENT(IN) :: nfill
      INTEGER, DIMENSION(nlnpe,nlpe,nlev),INTENT(IN) :: nnn
      INTEGER, DIMENSION(nlnpe,nlpe),     INTENT(IN) :: nbl
      INTEGER, DIMENSION(3,nclprc),       INTENT(IN) :: ipos

      REAL, DIMENSION(nc,nlnpe,nlpe),     INTENT(IN) :: emiss
      REAL, DIMENSION(nc,nlnpe,nlpe),     INTENT(IN) :: estore
      REAL, DIMENSION(nlnpe,nlpe,nlev),   INTENT(IN) :: so2em
      REAL, DIMENSION(nlnpe,nlpe,nlev),   INTENT(IN) :: be7em
      REAL, DIMENSION(nlnpe,nlpe,nlev),   INTENT(IN) :: be10em
      REAL, DIMENSION(nlnpe,nlpe,nlev),   INTENT(IN) :: lnoxem
      REAL, DIMENSION(nlnpe,nlpe,nlev),   INTENT(IN) :: acnoxem
      REAL, DIMENSION(nlonpe,nlatpe),     INTENT(IN) :: bl
      REAL, DIMENSION(nlonpe,nlatpe),     INTENT(IN) :: orog
      REAL, DIMENSION(4,nclprc),          INTENT(IN) :: pos
      REAL,                               INTENT(IN) :: time
      REAL,                               INTENT(IN) :: missings
      REAL,                               INTENT(IN) :: missingbe7
      REAL,                               INTENT(IN) :: missingbe10
      REAL,                               INTENT(IN) :: missinglnox
      REAL,                               INTENT(IN) :: missingacnox
      REAL, DIMENSION(nc+2,nclprc),       INTENT(OUT):: em

      INTEGER :: i
      INTEGER :: j
      INTEGER :: jj
      INTEGER :: k
      INTEGER :: l
      REAL :: ti
      REAL :: h
      REAL :: thet
      REAL,DIMENSION(nfill) :: zorog
      REAL,DIMENSION(nfill) :: zbl
      REAL :: rlong
      REAL :: rlat

      DO jj=1,nfill
! Find cells in boundary layer to recieve emissions
        i=ipos(1,jj)-lndat+1
        l=ipos(2,jj)-ltdat+1
        k=ipos(3,jj)
! DEPENDS ON: getmetpoint
        zorog(jj)=GETMETPOINT(pos(:,jj),orog,.TRUE.)
! DEPENDS ON: getmetpoint
        zbl(jj)=  GETMETPOINT(pos(:,jj),bl,.TRUE.)
        EM(I_LNOx,JJ)=0.
        EM(I_AcNOx,JJ)=0.
      END DO
      DO j=1,nc
        DO jj=1,nfill
! Find cells in boundary layer to recieve emissions
          i=ipos(1,jj)-lndat+1
          l=ipos(2,jj)-ltdat+1
          k=ipos(3,jj)
          em(j,jj)=0.
          IF (pos(4,jj)-earth_radius-zorog(jj) <= zbl(jj)               &
     &        .AND.j/=i_o3.AND.j/=i_hno3) THEN
            IF(nbl(i,l)<=0) THEN
!kk              WRITE(6,*) ' *** EGET: NBL=0 !',' at I=',i,' L=',l
!kk              WRITE(6,*) ' POS(4) = ',pos(4,jj)-Earth_Radius-Zorog(jj
!kk     &                   ' BL = ',Zbl(jj)
            ENDIF
            em(j,jj)=((emiss(j,i,l)+estore(j,i,l))/nbl(i,l))
          ENDIF
        END DO
      END DO

!kk   Computation of ti is loop invariant
! DEPENDS ON: secs
      ti=SECS(15.0,month,1)+MOD(time,daysec)

      DO jj=1,nfill
! Find cells in boundary layer to recieve emissions
        i=ipos(1,jj)-lndat+1
        l=ipos(2,jj)-ltdat+1
        k=ipos(3,jj)
! Set diurnal variation of isoprene emissions.
        IF (pos(4,jj)-earth_radius-zorog(jj) > zbl(jj)) THEN  !Not in BL
          em(i_c5h8,jj)=0.            ! Set isoprene emissions to zero
        ELSE                          ! In BL
! Get zenith angle for this time of day on 15th of the month
          rlat=(lat(l+ltdat-1-1)+lat(l+ltdat-1))/2.-90.
          IF(i+lndat-1<nlong) THEN
            rlong=(long(i+lndat-1)+long(i+lndat-1+1))/2.
          ELSE
            rlong=(long(i+lndat-1)+long(MOD(i+lndat-1,nlong)+1)+360.)/2.
          ENDIF
! DEPENDS ON: zen
          thet=ZEN(ti,rlat,rlong)
          IF(COS(thet)<0) THEN               ! It's dark
! emit any stored isoprene emissions
            em(i_c5h8,jj)=estore(i_c5h8,i,l)/nbl(i,l)
          ELSE                               ! It's light
! NB the stored isoprene emissions have already been multiplied
! by the appropriate cos(theta) for when they should have been
! emitted.
            em(i_c5h8,jj)=((emiss(i_c5h8,i,l)*cos(thet)+                &
     &                estore(i_c5h8,i,l))/nbl(i,l))
          ENDIF
        ENDIF

! Add in Stratospheric 'NO2', Lightning & Aircraft NOx,
!  SO2 volcanic, and Be emissions.  EM(I_LNOx) and EM(I_AcNOx) are
!  provided for flux accounting only.
        EM(I_NO,JJ)=EM(I_NO,JJ)+(MissingLNOx*LNOxEM(I,L,K)/NNN(I,L,K))
        EM(I_LNOx,JJ)=EM(I_LNOx,JJ)+(MissingLNOx*LNOxEM(I,L,K)/         &
     &                NNN(I,L,K))
        EM(I_NO,JJ)=EM(I_NO,JJ)+(MissingAcNOx*AcNOxEM(I,L,K)/NNN(I,L,K))
        EM(I_AcNOx,JJ)=EM(I_AcNOx,JJ)+(MissingAcNOx*AcNOxEM(I,L,K)/     &
     &                NNN(I,L,K))
        em(i_so2,jj)=em(i_so2,jj)+(missings*so2em(i,l,k)/nnn(i,l,k))
        em(i_be7,jj)= em(i_be7,jj)+                                     &
     &    (missingbe7* be7em(i,l,k) /nnn(i,l,k))
        em(i_be10,jj)=em(i_be10,jj)+                                    &
     &    (missingbe10*be10em(i,l,k)/nnn(i,l,k))
      END DO

      END SUBROUTINE EGET
#endif
