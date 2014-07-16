#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE LIGHTREAD(acnoxem,so2em,be7em,be10em,month,year1,      &
     &  totacnoxem,totso2em,totbe7em,totbe10em,p,lnp)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Read in aircraft, volcanic so2, and Be
!-                         emissions data
!-
!-   Inputs  : month, year, lbdat, p
!-   Outputs : acnoxem,totacnoxem,so2em,totso2em,
!-               be7em,totbe7em,be10em,totbe10em
!-   Controls:
!-
!
! Current Code Owner: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.5    07/10/94  Created.  D.S. Stevenson.
!  5.5    11/08/03  Fixed aircraft emission bug, still calculated
!                   1976/1984 emissions even if year was outside range.
!                   M.G. Sanderson
!  5.5    02/08/04  Added check for l_emiss_current. M.G. Sanderson
!  6.1    22/10/04  Function HEIGHT renamed ST_HEIGHT. M.G. Sanderson
!
!-
!VVV  v5.2.1 LIGHTREAD 17/VIII/01
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------

      INTEGER, INTENT(IN)               :: month
      INTEGER, INTENT(IN)               :: year1

      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: p
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: lnp
      REAL, DIMENSION(nlnpe,nlpe,nlev),        INTENT(OUT) :: acnoxem
      REAL, DIMENSION(nlnpe,nlpe,nlev),        INTENT(OUT) :: so2em
      REAL, DIMENSION(nlnpe,nlpe,nlev),        INTENT(OUT) :: be7em
      REAL, DIMENSION(nlnpe,nlpe,nlev),        INTENT(OUT) :: be10em
      REAL,                                    INTENT(OUT) :: totacnoxem
      REAL,                                    INTENT(OUT) :: totso2em
      REAL,                                    INTENT(OUT) :: totbe7em
      REAL,                                    INTENT(OUT) :: totbe10em

      INTEGER :: j
      INTEGER :: jj
      INTEGER :: k
      INTEGER :: kk
      INTEGER :: kktop
      INTEGER :: kkbot
      INTEGER :: year
      INTEGER :: month1
      INTEGER :: month2

      REAL :: f
      REAL :: f1
      REAL :: etatop
      REAL :: etabot
      REAL :: area
      REAL :: p1
      REAL :: ptop
      REAL :: pbot
      REAL, DIMENSION(4) :: rpos
      REAL, DIMENSION(nlong,mnlat,20)   :: aircraft
      REAL, DIMENSION(nlong,mnlat,20)   :: aircraft_1976
      REAL, DIMENSION(nlong,mnlat,20)   :: aircraft_1984
      REAL, DIMENSION(nlong,mnlat,20)   :: aircraft_1992
      REAL, DIMENSION(nlong,mnlat,20)   :: aircraft_2015
      REAL, DIMENSION(nlong,mnlat,20)   :: aircraft_2050
      REAL, DIMENSION(nlong,mnlat,20)   :: aircraft_one
      REAL, DIMENSION(nlong,mnlat,20)   :: aircraft_two
      REAL, DIMENSION(mnlat,34)         :: be7
      REAL, DIMENSION(mnlat,34)         :: be10
      REAL, DIMENSION(9,34)             :: be7in
      REAL, DIMENSION(9,34)             :: be10in
      REAL, DIMENSION(35)               :: depth
      REAL, DIMENSION(35)               :: press
      REAL, DIMENSION(nlong,mnlat)      :: so2
      REAL, DIMENSION(nlong,mnlat,nlev) :: sumno2
      REAL, PARAMETER :: p_ref=1.0e5         ! Ref. Surface Press

      CHARACTER(LEN=2)                  :: smonth
      CHARACTER(LEN=2)                  :: smonth1
      CHARACTER(LEN=2)                  :: smonth2
      CHARACTER(LEN=1)                  :: slat

      COMMON /LIGHT/ aircraft

! Need to set all the variables to zero?
      aircraft_1976 = 0.0
      aircraft_1984 = 0.0
      aircraft_1992 = 0.0
      aircraft_2015 = 0.0
      aircraft_2050 = 0.0
      sumno2 = 0.0

      IF (l_emiss_current) THEN
        year = year1
      ELSE
        year = emiss_year
        WRITE(6,*) 'Warning: Aircraft emissions are for ',emiss_year
      END IF

! N.B. NO2 emissions are in kg[NO2]/day, not kg[N]
!
! NASA data from HSRP final version (V971020), interpolated onto
! STOCHEM grid.
! NASA 1992 and 2015 data is monthly, 1976 and 1984 are for Feb, May,
! August and November,  2050 data is an annual average.

! Write month digit(s) into string
      WRITE(smonth,'(I2.2)') month

      IF (aircrafton) THEN
! Aircraft NO2 data, kg/day/grid square
! Monthly:
        IF (year > 1984 .AND. year <= 2015) THEN
          OPEN(21,FILE=TRIM(emdir2)//'dt092'//smonth//'k_9672.v97',     &
     &      STATUS='OLD')
          READ(21,'(6E13.5)') aircraft_1992(:,:,1:20)
          CLOSE(21)
        END IF

        IF (year > 1992 .AND. year <= 2050) THEN
          OPEN(21,FILE=TRIM(emdir2)//'dt015'//smonth//'k_9672.v97',     &
     &      STATUS='OLD')
          READ(21,'(6E13.5)') aircraft_2015(:,:,1:20)
          CLOSE(21)
        END IF

! Data available 4 times per year for 1976 and 1984
        IF (month == 2 .OR. month == 5 .OR. month == 8 .OR. month == 11)&
     &    THEN
! Read one month
          IF (year > 1950 .AND. year <= 1984) THEN
            OPEN(21,FILE=TRIM(emdir2)//'dt076'//smonth//'k_9672.v97',   &
     &        STATUS='OLD')
            READ(21,'(6e13.5)') aircraft_1976(:,:,1:20)
            CLOSE(21)
          END IF

          IF (year > 1976 .AND. year <= 1992) THEN
            OPEN(21,FILE=TRIM(emdir2)//'dt084'//smonth//'k_9672.v97',   &
     &        STATUS='OLD')
            READ(21,'(6e13.5)') aircraft_1984(:,:,1:20)
            CLOSE(21)
          END IF
        ELSE
! Read two months, then interpolate
          IF (month == 1 .or. month == 4 .or. month == 7 .OR.           &
     &      month == 10) then
            month1 = month - 2
            month2 = month + 1
          ELSE IF (month == 3 .or. month == 6 .or. month == 9 .OR.      &
     &      month == 12) THEN
            month1=month-1
            month2=month+2
          END IF

          IF (month1 < 1) THEN
            WRITE(smonth1,'(i2.2)') month1+12
          ELSE
            WRITE(smonth1,'(i2.2)') month1
          END IF
          IF (month2 > 12) THEN
            WRITE(smonth2,'(i2.2)') month2-12
          ELSE
            WRITE(smonth2,'(i2.2)') month2
          END IF

          IF (year > 1950 .AND. year <= 1984) THEN
            OPEN(21,FILE=TRIM(emdir2)//'dt076'//smonth1//'k_9672.v97',  &
     &        STATUS='OLD')
            READ(21,'(6e13.5)') aircraft_one(:,:,1:20)
            CLOSE(21)

            OPEN(21,FILE=TRIM(emdir2)//'dt076'//smonth2//'k_9672.v97',  &
     &        STATUS='OLD')
            READ(21,'(6e13.5)') aircraft_two(:,:,1:20)
            CLOSE(21)
            aircraft_1976 = aircraft_one +                              &
     &        ((aircraft_two-aircraft_one)/3.0) * REAL(month-month1)
          END IF

          IF (year > 1976 .AND. year <= 1992) THEN
            OPEN(21,FILE=TRIM(emdir2)//'dt084'//smonth1//'k_9672.v97',  &
     &        STATUS='OLD')
            READ(21,'(6e13.5)') aircraft_one(:,:,1:20)
            CLOSE(21)

            OPEN(21,FILE=TRIM(emdir2)//'dt084'//smonth2//'k_9672.v97',  &
     &        STATUS='OLD')
            READ(21,'(6e13.5)') aircraft_two(:,:,1:20)
            CLOSE(21)
              aircraft_1984 = aircraft_one +                            &
     &          ((aircraft_two-aircraft_one)/3.0) * REAL(month-month1)
          END IF
        END IF

! Annual data, 2050 scenario IS92a with fuel efficiency emphasis
        IF (year > 2015) THEN
          OPEN(21,FILE=TRIM(emdir2)//'dt050I0k_9672.v97',STATUS='OLD')
          READ(21,'(6e13.5)') aircraft_2050(:,:,1:20)
          CLOSE(21)
        END IF
!
      END IF   ! IF (aircrafton)

! Interpolate aircraft for year (zero b4 1950; 2050 values used after 20
      IF (.NOT. aircrafton .OR. scenario == 'pi') THEN
        aircraft = 0.0         ! zero aircraft
      ELSE
        IF (year <= 1950) THEN
          aircraft = 0.0
        ELSE IF (year > 1950 .AND. year <= 1976) THEN
          f = REAL(1966-year) / (1976.0-1950.0)
          aircraft=(1.0-f)*aircraft_1976
        ELSE IF (year > 1976 .AND. year <= 1984) THEN
          f = REAL(1984-year) / (1984.0-1976.0)
          aircraft = f*aircraft_1976 + (1.0-f)*aircraft_1984
        ELSE IF (year > 1984 .AND. year <= 1992) THEN
          f = REAL(1992-year) / (1992.0-1984.0)
          aircraft = f*aircraft_1984 + (1.0-f)*aircraft_1992
        ELSE IF (year > 1992.and.year <= 2015) THEN
          f = REAL(2015-year) / (2015.0-1992.0)
          aircraft = f*aircraft_1992 + (1.0-f)*aircraft_2015
        ELSE IF (year > 2015.and.year <= 2050) THEN
          f = REAL(2050-year) / (2050.0-2015.0)
          aircraft = f*aircraft_2015 + (1.0-f)*aircraft_2050
        ELSE IF(year > 2050) THEN
          aircraft = aircraft_2050
        END IF
      END IF

      WRITE(6,*) 'Aircraft emissions: ',sum(aircraft),' kg[NO2]/day'

! Convert from per day to per second
      SUMNO2 = aircraft / daysec

! Volcanic SO2 data, kg[S]/yr/grid square (2-D).
! Continous
      OPEN(21,FILE=TRIM(emdir2)//'so2volc_cont_96_72.dat',STATUS='OLD')
      READ(21,'(6e13.5)') so2
      CLOSE(21)
! Sporadic
!!      OPEN(21,FILE=TRIM(emdir2)//'so2volc_sporad_96_72.dat',STATUS='OL
!!      READ(21,'(6e13.5)') SO2
!!      CLOSE(21)
!   Distribute in vertical to 300mb = level 14
! DEPENDS ON: st_height
      kk = ST_HEIGHT(0.3,'Eta_stochem')       ! 0.3 ~200mb
      so2em = 0.0
      DO k=1,kk
        so2em(:,:,k) = so2(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1) /    &
     &    (REAL(kk)*ysec)
      END DO

! Read in Be production in atoms/g(air)
! Read in column depths in g(air)/cm2
      DO j=1,9
        WRITE(slat,'(i1)') j
        OPEN(21,FILE=TRIM(emdir2)//'be7'//slat//'0550.out',STATUS='OLD')
        OPEN(22,FILE=TRIM(emdir2)//'be'//slat//'0550.out',STATUS='OLD')
        DO k=1,34
          READ(21,*) depth(k), be7in(j,k)
          READ(22,*) depth(k), be10in(j,k)
        END DO
        CLOSE(21)
        CLOSE(22)
      END DO

      depth(35) = p_ref * 0.1 / g  ! Set bottom depth to P_REF in g/cm2
      press = depth*10.0*g         ! convert from g/cm2 to N/m2

      DO jj=1,mnlat
! Calculate J index from latitude at centre of grid square
        j = INT(ABS((lat(jj-1)+lat(jj))/2.0 - 90.0)/10) + 1
        area = (2.0*pi*earth_radius**2) *                               &
     &    (SIN((lat(jj)-90.)*pi_over_180) -                             &
     &    SIN((lat(jj-1)-90.)*pi_over_180)) / REAL(nlong)

! Multiply by layer depth and area to go from g/g(air) to g/grid volume
        be7(jj,:) = be7in(j,:) *1.0e4*(depth(2:35)-depth(1:34))*area
        be10(jj,:)= be10in(j,:)*1.0e4*(depth(2:35)-depth(1:34))*area
      END DO

! Interpolation onto height grid
      be7em = 0.0
      be10em = 0.0

! to get pressure, take a point on Greenwich Meridien in the middle of t
! latitude range for the processor
      rpos(1) = (longm_half(lnbound)+longm_half(lnbound+nlonpe-1))/2.0
      rpos(2) = (latm_half(lobound)+latm_half(lobound+nlatpe-2))/2.0
      DO k=1,34
        IF (k == 1) THEN
          etatop = eta_stochem(nlev)
!         rpos(1)=(longm_half(lnbound)+longm_half(lnbound+nlonpe-1))/2.0
!         rpos(2)=(latm_half(lobound)+latm_half(lobound+nlatpe-2))/2.0
          rpos(3) = etatop
! DEPENDS ON: eta2p
          ptop=ETA2P(rpos,lnp)
        ELSE
!         rpos(1)=(longm_half(lnbound)+longm_half(lnbound+nlonpe-1))/2.0
!         rpos(2)=(latm_half(lobound)+latm_half(lobound+nlatpe-2))/2.0
          rpos(3) = 0.0
! DEPENDS ON: getmetpoint
          IF (GETMETPOINT(rpos,p(:,:,0),.TRUE.) > press(k)) THEN
! DEPENDS ON: p2eta
            etatop = P2ETA(press(k),rpos,p)
          ELSE
            etatop = 0.0
          END IF
          ptop = press(k)
        END IF
        IF (k == 34) THEN
          etabot = 0.0
          pbot = p_ref
        ELSE
!         rpos(1)=(longm_half(lnbound)+longm_half(lnbound+nlonpe-1))/2.0
!         rpos(2)=(latm_half(lobound)+latm_half(lobound+nlatpe-2))/2.0
          rpos(3) = 0.0
! DEPENDS ON: getmetpoint
          IF (GETMETPOINT(rpos,p(:,:,0),.TRUE.) > press(k+1)) THEN
! DEPENDS ON: p2eta
            etabot = P2ETA(press(k+1),rpos,p)
          ELSE
            etabot = 0.0
          END IF
          pbot = press(k+1)
        END IF
! DEPENDS ON: st_height
        kktop = ST_HEIGHT(etatop,'Eta_stochem')
! DEPENDS ON: st_height
        kkbot = ST_HEIGHT(etabot,'Eta_stochem')
!       rpos(1)=(longm_half(lnbound)+longm_half(lnbound+nlonpe-1))/2.0
!       rpos(2)=(latm_half(lobound)+latm_half(lobound+nlatpe-2))/2.0
        rpos(3) = eta_stochem(kktop)
! DEPENDS ON: eta2p
        p1 = ETA2P(rpos,lnp)

! F1 is the fraction of emission layer K lying within grid layer KKTOP
        f1 = MIN((p1-ptop)/(pbot-ptop),1.0)
        be7em(:,:,kktop) = be7em(:,:,kktop) +                           &
     &    SPREAD(be7(ltdat:ltdat+nlpe-1,k),1,nlnpe) * f1
        be7em(:,:,kkbot) = be7em(:,:,kkbot) +                           &
     &    SPREAD(be7(ltdat:ltdat+nlpe-1,k),1,nlnpe) * (1.0-f1)
        be10em(:,:,kktop) = be10em(:,:,kktop) +                         &
     &    SPREAD(be10(ltdat:ltdat+nlpe-1,k),1,nlnpe) * f1
        be10em(:,:,kkbot) = be10em(:,:,kkbot) +                         &
     &    SPREAD(be10(ltdat:ltdat+nlpe-1,k),1,nlnpe) * (1.0-f1)
      END DO

!   Sum and convert to (       molecules / grid square       )
!                      ( ----------------------------------- )  /s
!                      ( molecules air / lagrangian particle )

      totbe7em = SUM(be7)*nlong / lmolec
      totbe10em = SUM(be10)*nlong / lmolec
      be7em = be7em / lmolec
      be10em = be10em / lmolec
      acnoxem = sumno2(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1,:)*na /   &
     &  (mno2*lmolec)
      totacnoxem = SUM(sumno2)*na / (mno2*lmolec)
      so2em = so2em*na / (msul*lmolec)
      so2 = so2*na / (msul*lmolec)
      totso2em = SUM(so2)

      END SUBROUTINE LIGHTREAD
#endif
