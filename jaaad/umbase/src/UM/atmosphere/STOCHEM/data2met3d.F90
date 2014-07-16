#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE DATA2MET3D(x,y,mw)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Convert field from 3-D STOCHEM grid to
!-                         3-D on met Eta_theta grid by
!-                         linear interpolation, and convert from
!-                         vmr to mmr.
!-
!-   Inputs  : X
!-   Outputs : Y
!-   Controls:
!-
!
! Current Code Owner: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  5.5    07/02/03  Created.  C.E. Johnson
!  5.5    15/05/03  Now creates a 3D field. W.J. Collins
!  6.1    22/10/04  No change.
!
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!-----------------------------------------------------------------------

      REAL, DIMENSION(nlnpe,nlpe,nlev),        INTENT(IN)  :: x
      REAL,                                    INTENT(IN)  :: mw
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(OUT)  :: y

      INTEGER :: i
      INTEGER :: i2
      INTEGER :: j
      INTEGER :: jj
      INTEGER :: k
      INTEGER :: l

      REAL :: Delta
      REAL :: Delta2
      REAL, DIMENSION(nlong,mnlat,nlev)     :: global
      REAL, DIMENSION(nlong,mnlat,nmetlev)  :: intrp1
      REAL, DIMENSION(nlong,nlatpe,nmetlev) :: intrp2
      REAL, DIMENSION(nlong,nmetlev)        :: row0
      REAL, DIMENSION(nlong,nmetlev)        :: temporary
      REAL, DIMENSION(nlev)                 :: Eta_mid
      REAL, DIMENSION(nlong)                :: Lon_mid
      REAL, DIMENSION(mnlat)                :: Lat_mid
      REAL                                  :: lon_h
      REAL                                  :: lat_h

      LOGICAL                               :: l_polar

      COMMON /SUM3_D/ global

! DEPENDS ON: sum3d
      CALL SUM3D(x)

! Interpolate from mid point of STOCHEM grids,
!  so make appropriate ones
      Eta_mid=(Eta_Stochem(1:nlev)+Eta_Stochem(0:nlev-1))/2.0
      Lat_mid=(lat(1:mnlat)+lat(0:mnlat-1))/2.0
      Lon_mid(1:nlong-1)=(long(1:nlong-1)+long(2:nlong))/2.0
      lon_mid(nlong)=(long(nlong)+long(1)+360.)/2.0

! Interpolate with height first
      DO k=1,nmetlev
        IF(Eta_theta(k) < Eta_mid(1)) THEN
          intrp1(:,:,k)=Global(:,:,1)
        ELSEIF(Eta_theta(k) > Eta_mid(nlev)) THEN
          intrp1(:,:,k)=Global(:,:,nlev)
        ELSE
          l=BISECT(eta_mid,Eta_theta(k),nlev)
          Delta=Eta_mid(l+1)-Eta_mid(l)
          Delta2=Eta_Theta(k)-Eta_mid(l)
          intrp1(:,:,k)=Global(:,:,l)+(Global(:,:,l+1)-                 &
     &                Global(:,:,l))*Delta2/Delta
        ENDIF
      ENDDO

      DO j=1,nlatpe
        jj=j+lobound-1
        lat_h=Latm_Half(jj)
        l_polar=.FALSE.
        IF(lat_h>180.) THEN
          lat_h=360.-lat_h
          l_polar=.TRUE.
        ELSEIF(lat_h<0.) THEN
          lat_h=-lat_h
          l_polar=.TRUE.
        ENDIF
        Delta=dlat
        IF(Lat_H < Lat_mid(1)) THEN
          Delta2=Lat_H-Lat_mid(1)+dlat
! row0 is row 1 of intrp1 shifted by 180 degrees
          row0(1:nlong/2,:)=intrp1(nlong/2+1:nlong,1,:)
          row0(nlong/2+1:nlong,:)=intrp1(1:nlong/2,1,:)
          intrp2(:,j,:)=row0+(intrp1(:,1,:)-row0)*Delta2/Delta
        ELSEIF(Lat_H > Lat_mid(mnlat)) THEN
          Delta2=Lat_H-Lat_mid(mnlat)
! row0 is row mnlat of intrp1 shifted by 180 degrees
          row0(1:nlong/2,:)=intrp1(nlong/2+1:nlong,mnlat,:)
          row0(nlong/2+1:nlong,:)=intrp1(1:nlong/2,mnlat,:)
          intrp2(:,j,:)=intrp1(:,mnlat,:)+(row0                         &
     &                  -intrp1(:,mnlat,:))*Delta2/Delta
        ELSE
          i=BISECT(Lat_mid,Lat_H,mnlat)
          Delta2=Lat_H-Lat_mid(i)
          intrp2(:,j,:)=intrp1(:,i,:)+(intrp1(:,i+1,:)-intrp1(:,i,:))   &
     &             *Delta2/Delta
        ENDIF
        IF(l_polar) THEN
! if we are dealing with a halo beyond the poles, shift row by 180
          temporary(1:nlong/2,:)=intrp2(nlong/2+1:nlong,j,:)
          temporary(nlong/2+1:nlong,:)=intrp2(1:nlong/2,j,:)
          intrp2(:,j,:)=temporary
        ENDIF
      ENDDO

      Delta=dlong
      DO j=1,nlonpe
        jj=j+lnbound-1
        lon_h=AMOD(longm_half(jj)+360.,360.)
        IF(lon_h < Lon_mid(1)) THEN
          i=nlong
          i2=1
          delta2=lon_h+360.-lon_mid(nlong)
          y(j,:,:)=intrp2(i,:,:)+(intrp2(i2,:,:)-intrp2(i,:,:))         &
     &           *Delta2/Delta
        ELSEIF(lon_h > Lon_mid(nlong)) THEN
          i=nlong
          i2=1
          delta2=lon_h-lon_mid(nlong)
          y(j,:,:)=intrp2(i,:,:)+(intrp2(i2,:,:)-intrp2(i,:,:))         &
     &           *Delta2/Delta
        ELSE
          i=BISECT(Lon_mid,lon_h,nlong)
          i2=i+1
          delta2=lon_h-Lon_mid(i)
          y(j,:,:)=intrp2(i,:,:)+(intrp2(i2,:,:)-intrp2(i,:,:))         &
     &           *Delta2/Delta
        ENDIF
      ENDDO

! Convert into mass mixing ratio for the UM:
      y(:,:,:)=y(:,:,:)*mw/mair

      CONTAINS
      INTEGER FUNCTION BISECT(earray,z,ji)
! Bisection method

      REAL, DIMENSION(*), INTENT(IN) :: earray
      REAL,    INTENT(IN) :: z
      INTEGER, INTENT(IN) :: ji

      INTEGER :: jm
      INTEGER :: jl
      INTEGER :: ju

      jl=1
      ju=ji
      jm=(ju+jl)/2
      DO
        IF (ju-jl == 1) EXIT
        IF (z > earray(jm)) THEN
          jl=jm
        ELSE
          ju=jm
        ENDIF
        jm=(ju+jl)/2
      ENDDO
      BISECT=jm

      END FUNCTION BISECT

      END SUBROUTINE DATA2MET3D
#endif
