#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EMREAD(emiss,year,class)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : READ IN EMISSIONS DATA
!-    reads in emissions in tonnes/year per grid square
!-    converts to molecules/s per grid square (output EMISS array)
!-    non-seasonal components
!-
!-   Inputs  : CLASS,YEAR
!-   Outputs : EMISS
!-   Controls:
!-
!
! Current Code Owner: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.5    15/12/93  Created.  W.J. Collins
!  5.1    11/12/00  Units now in kg not g. C.E. Johnson
!  5.5    13/02/04  EMDIR redefined so need TRIM function. K. Ketelsen
!  6.1    06/09/04  Now selects year using l_emiss_current.
!                   M.G. Sanderson
!  6.2    01/03/06  Corrected bug in calculation of correct year
!                   for emissions.  C.E. Johnson.
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER,                        INTENT(IN)  :: year
      REAL, DIMENSION(5,nc),          INTENT(IN)  :: class
      REAL, DIMENSION(nc,nlnpe,nlpe), INTENT(OUT) :: emiss

      INTEGER :: iy1
      INTEGER :: iy2
      INTEGER :: inc
      INTEGER :: yearx

      REAL :: rinc
      REAL :: s_anth_nmhc
      REAL, DIMENSION(nlong,mnlat) :: anth_nox_year1
      REAL, DIMENSION(nlong,mnlat) :: anth_nox_year2
      REAL, DIMENSION(nlong,mnlat) :: anth_nox
      REAL, DIMENSION(nlong,mnlat) :: anth_nox_ship
      REAL, DIMENSION(nlong,mnlat) :: anth_co_year1
      REAL, DIMENSION(nlong,mnlat) :: anth_co_year2
      REAL, DIMENSION(nlong,mnlat) :: anth_co
      REAL, DIMENSION(nlong,mnlat) :: anth_ch4_year1
      REAL, DIMENSION(nlong,mnlat) :: anth_ch4_year2
      REAL, DIMENSION(nlong,mnlat) :: anth_ch4
      REAL, DIMENSION(nlong,mnlat) :: anth_nmhc_year1
      REAL, DIMENSION(nlong,mnlat) :: anth_nmhc_year2
      REAL, DIMENSION(nlong,mnlat) :: anth_nmhc
      REAL, DIMENSION(nlong,mnlat) :: anth_so2_year1
      REAL, DIMENSION(nlong,mnlat) :: anth_so2_year2
      REAL, DIMENSION(nlong,mnlat) :: anth_so2
      REAL, DIMENSION(nlong,mnlat) :: anth_nh3

      CHARACTER(LEN=4) :: cy1
      CHARACTER(LEN=4) :: cy2
      CHARACTER(LEN=7) :: fname1
      CHARACTER(LEN=72):: cmessage

      IF (scenario /= 'pi') THEN    ! anthropogenic components included

! Select the scenario, INC is the years between sucessive files.

! Find datasets to form interpolation pair (SRES).
        inc = 10
        IF (.NOT. l_emiss_current) THEN
          yearx = emiss_year
          WRITE(6,*) 'EMREAD: Using scenario emissions from: ',         &
     &      emiss_year
        ELSE
          yearx = year
        END IF
        iy1 = (yearx / 10) * 10
        IF (iy1 < 1990 .OR. iy1 > 2090) THEN    ! Stop if out of range
          cmessage = 'Emission year out of range for scenario'
          WRITE(6,*) 'EMREAD: ', cmessage, iy1
! DEPENDS ON: ereport
          CALL EREPORT('Emread',1,cmessage)
        END IF
        iy2 = iy1 + inc
        rinc = REAL(inc)

! For IIASA CLE scenario, only have years 2000 and 2030. Set iy2 to iy1
! this scenario selected
        IF (scenario == 'bu' .OR. scenario == 'mf' .OR.                 &
     &    scenario == 'A2') iy2 = iy1

! Write into character variables.
        WRITE(cy1,'(i4)') iy1
        WRITE(cy2,'(i4)') iy2

! NOx data, kg(N)/yr per grid square
        OPEN(21,FILE=trim(emdir)//scenario//cy1//'NOx_96_72.dat',       &
     &    STATUS='OLD')
        READ(21,*) anth_nox_year1
        CLOSE(21)
        OPEN(21,FILE=trim(emdir)//scenario//cy2//'NOx_96_72.dat',       &
     &    STATUS='OLD')
        READ(21,*) anth_nox_year2
        CLOSE(21)
        anth_nox = anth_nox_year1 + real(yearx-iy1) * (anth_nox_year2 - &
     &    anth_nox_year1) / rinc

! Read marine emissions kg (N)/yr and add to anthropogenic
        OPEN(21,FILE=trim(emdir)//'ship1990_NOx_96_72.dat',             &
     &    STATUS='OLD')
        READ(21,*) anth_nox_ship
        CLOSE(21)

! CO  data, kg(CO)/yr per grid square
        OPEN(21,FILE=trim(emdir)//scenario//cy1//'CO_96_72.dat',        &
     &    STATUS='OLD')
        READ(21,*) anth_co_year1
        CLOSE(21)
        OPEN(21,FILE=trim(emdir)//scenario//cy2//'CO_96_72.dat',        &
     &    STATUS='OLD')
        READ(21,*) anth_co_year2
        CLOSE(21)
        anth_co = anth_co_year1 + REAL(yearx-iy1) * (anth_co_year2 -    &
     &    anth_co_year1) / rinc

! CH4 data, kg(CH4)/yr per grid square
! (includes animals & paddies)
        OPEN(21,FILE=trim(emdir)//scenario//cy1//'CH4_96_72.dat',       &
     &    STATUS='OLD')
        READ(21,*) anth_ch4_year1
        CLOSE(21)
        OPEN(21,FILE=trim(emdir)//scenario//cy2//'CH4_96_72.dat',       &
     &    STATUS='OLD')
        READ(21,*) anth_ch4_year2
        CLOSE(21)
        anth_ch4 = anth_ch4_year1 + REAL(yearx-iy1) * (anth_ch4_year2 - &
     &    anth_ch4_year1) / rinc

! NMHC data, kg(NMHC)/yr per grid square
        OPEN(21,FILE=trim(emdir)//scenario//cy1//'NMVOC_96_72.dat',     &
     &    STATUS='OLD')
        READ(21,*) anth_nmhc_year1
        CLOSE(21)
        OPEN(21,FILE=trim(emdir)//scenario//cy2//'NMVOC_96_72.dat',     &
     &    STATUS='OLD')
        READ(21,*) anth_nmhc_year2
        CLOSE(21)
        anth_nmhc = anth_nmhc_year1 + REAL(yearx-iy1) *                 &
     &    (anth_nmhc_year2 - anth_nmhc_year1) / rinc

! SO2 data, kg(S)/yr per grid square
        OPEN(21,FILE=trim(emdir)//scenario//cy1//'SO2_96_72.dat',       &
     &    STATUS='OLD')
        READ(21,*) anth_so2_year1
        CLOSE(21)
        OPEN(21,FILE=trim(emdir)//scenario//cy2//'SO2_96_72.dat',       &
     &    STATUS='OLD')
        READ(21,*) anth_so2_year2
        CLOSE(21)
        anth_so2 = anth_so2_year1 + REAL(yearx-iy1) * (anth_so2_year2 - &
     &    anth_so2_year1) / rinc

! NH3 data, kg(NH3)/yr per grid square
        OPEN(21,FILE=trim(emdir)//'nh3ann_96_72.dat',STATUS='OLD')
        READ(21,*) anth_nh3
        CLOSE(21)

! Calculate anthropogenic emissions
        s_anth_nmhc = SUM(anth_nmhc)
        emiss(i_no,:,:) =((class(1,i_no)/SUM(anth_nox))*                &
     &    anth_nox(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+             &
     &    (class(5,I_NO)/SUM(anth_nox_ship))*                           &
     &    anth_nox_ship(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*       &
     &    na/(mnit*ysec)
        emiss(i_co,:,:) =(class(1,i_co)/SUM(anth_co))*                  &
     &    anth_co(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*              &
     &    na/(mco*ysec)
        emiss(i_ch4,:,:)=(class(1,i_ch4)/SUM(anth_ch4))*                &
     &    anth_ch4(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*             &
     &    na/(mch4*ysec)
        emiss(i_hcho,:,:)=(class(1,i_hcho)/s_anth_nmhc)*                &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(mhcho*ysec)
        emiss(i_c2h6,:,:)=(class(1,i_c2h6)/s_anth_nmhc)*                &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(mc2h6*ysec)
        emiss(i_ch3cho,:,:)=(class(1,i_ch3cho)/s_anth_nmhc)*            &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(mch3cho*ysec)
        emiss(i_nc4h10,:,:)=(class(1,i_nc4h10)/s_anth_nmhc)*            &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(mnc4h10*ysec)
        emiss(i_so2,:,:)=(class(1,i_so2)/SUM(anth_so2))*                &
     &    anth_so2(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*             &
     &    na/(msul*ysec)
        emiss(i_c2h4,:,:)=(class(1,i_c2h4)/s_anth_nmhc)*                &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(mc2h4*ysec)
        emiss(i_c3h6,:,:)=(class(1,i_c3h6)/s_anth_nmhc)*                &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(mc3h6*ysec)
        emiss(i_c3h8,:,:)=(class(1,i_c3h8)/s_anth_nmhc)*                &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(mc3h8*ysec)
        emiss(i_ch3oh,:,:)=(class(1,i_ch3oh)/s_anth_nmhc)*              &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(mch3oh*ysec)
        emiss(i_acetone,:,:)=(class(1,i_acetone)/s_anth_nmhc)*          &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(macetone*ysec)
        emiss(i_oxyl,:,:)=(class(1,i_oxyl)/s_anth_nmhc)*                &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(moxyl*ysec)
        emiss(i_toluen,:,:)=(class(1,i_toluen)/s_anth_nmhc)*            &
     &    anth_nmhc(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*            &
     &    na/(mtoluen*ysec)
        emiss(i_nh3,:,:)=(class(1,i_nh3)/SUM(anth_nh3))*                &
     &    anth_nh3(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*             &
     &    na/(mnit*ysec)
! H2 (anth) using co distribution
        emiss(i_h2,:,:)=(class(1,i_h2)/SUM(anth_co))*                   &
     &    anth_co(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*              &
     &    na/(mh2*ysec)

      END IF       ! anthropogenic components

! Write out the CLASS array for key species

      WRITE(6,*) 'Emissions of key species (kg):'
      WRITE(6,*) 'Species  Anthrop  Bio-Burn  Vegetation  Soils   Ocean'
      WRITE(6,'(A6,5(1PE12.2))') ' NO   ',class(:,i_no)
      WRITE(6,'(A6,5(1PE12.2))') ' CO   ',class(:,i_co)
      WRITE(6,'(A6,5(1PE12.2))') ' CH4  ',class(:,i_ch4)
      WRITE(6,'(A6,5(1PE12.2))') ' H2   ',class(:,i_h2)
      WRITE(6,'(A6,5(1PE12.2))') ' SO2  ',class(:,i_so2)
      WRITE(6,'(A6,5(1PE12.2))') ' C2H4 ',class(:,i_c2h4)
      WRITE(6,'(A6,5(1PE12.2))') ' C2H6 ',class(:,i_c2h6)
      WRITE(6,'(A6,5(1PE12.2))') ' C3H6 ',class(:,i_c3h6)
      WRITE(6,'(A6,5(1PE12.2))') ' C3H8 ',class(:,i_c3h8)
      WRITE(6,'(A6,5(1PE12.2))') ' CH3OH',class(:,i_ch3oh)
      WRITE(6,'(A6,5(1PE12.2))') ' C5H8 ',class(:,i_c5h8)

      END SUBROUTINE EMREAD
#endif
