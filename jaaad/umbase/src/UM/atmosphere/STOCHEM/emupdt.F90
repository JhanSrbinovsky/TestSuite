#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EMUPDT(emiss,class,month,year)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : READ IN EMISSIONS DATA
!-
!-   Inputs  :
!-   Outputs : EMISS
!-   Controls:
!-
!
! Current Code Owner: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.5    10/01/94  Created.  W.J. Collins
!  4.5    23/10/98  Bio burn from GEIA (Cooke+Wilson 96) D.S. Stevenson
!  5.1    11/12/00  Units now in kg not g. C.E. Johnson
!  5.5    13/02/04  EMDIR redefined so need TRIM function. K. Ketelsen
!  6.1    22/10/04  Now reads in wetland emissions from Fung et al.
!                   JGR (1991), original files have anomaly over
!                   S.America. M.G. Sanderson
!  6.2    01/03/06  Reads in different biomass burning emission files
!                   for certain scenarios.  M.G. Sanderson.
!
!-
!VVV  v4.5   EMUPDT 11/12/00 - CLASS in kg, not g. Molecular masses
!                              specified as constants (kg/mol)
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD, soil_tile => soil
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER,                           INTENT(IN) :: month
      INTEGER,                           INTENT(IN) :: year
      REAL, DIMENSION(5,nc),             INTENT(IN) :: class
      REAL, DIMENSION(nc,nlnpe,nlpe), INTENT(INOUT) :: emiss

      REAL, DIMENSION(nlong,mnlat) :: burn
      REAL, DIMENSION(nlong,mnlat) :: gfed_nox
      REAL, DIMENSION(nlong,mnlat) :: gfed_ch4
      REAL, DIMENSION(nlong,mnlat) :: gfed_co
      REAL, DIMENSION(nlong,mnlat) :: gfed_nmhc
      REAL, DIMENSION(nlong,mnlat) :: gfed_so2
      REAL, DIMENSION(nlong,mnlat) :: gfed_nh3
      REAL, DIMENSION(nlong,mnlat) :: veg
      REAL, DIMENSION(nlong,mnlat) :: soil
      REAL, DIMENSION(nlong,mnlat) :: ocean
      REAL, DIMENSION(nlong,mnlat) :: wetland
!     REAL, DIMENSION(nlong,mnlat) :: isop
      REAL, DIMENSION(nlong,mnlat) :: ocean_dms
      REAL, DIMENSION(nlong,mnlat) :: radon
!     REAL, DIMENSION(251)         :: bb_factor
!     INTEGER :: bb_year

! Set to true if CLE, MFR or A2_  scenarios selected.
      LOGICAL :: l_bu_selected
      CHARACTER(LEN=2) :: smonth

! wetl_fac increases wetland emissions to 158 Tg/yr. The data files
! have units of kg grid cell-1 yr-1.
      REAL :: wetl_fac = 158.0 / 111.901

      l_bu_selected = (scenario == 'bu' .OR. scenario == 'mf' .OR.      &
     &  scenario == 'A2')

! Write month digit(s) into string and read monthly data
! All except biomass (normalised data) are in kg/grid cell.
      WRITE(smonth,'(I2.2)') MONTH

      IF (.NOT.l_bu_selected) THEN
        OPEN(21,FILE=trim(emdir)//'biomass'//smonth//'_96_72.dat',      &
     &    STATUS='OLD')
        READ(21,*) burn       ! normalised data
        CLOSE(21)
      ELSE
        OPEN(21,FILE=trim(emdir)//'bb_GFED_CH4_96_72.dat'//smonth,      &
     &    STATUS='OLD')
        READ(21,*) gfed_ch4       ! normalised data
        CLOSE(21)
        OPEN(21,FILE=trim(emdir)//'bb_GFED_CO_96_72.dat'//smonth,       &
     &    STATUS='OLD')
        READ(21,*) gfed_co       ! normalised data
        CLOSE(21)
        OPEN(21,FILE=trim(emdir)//'bb_GFED_NOX_96_72.dat'//smonth,      &
     &    STATUS='OLD')
        READ(21,*) gfed_nox       ! normalised data
        CLOSE(21)
        OPEN(21,FILE=trim(emdir)//'bb_GFED_NMHC_96_72.dat'//smonth,     &
     &    STATUS='OLD')
        READ(21,*) gfed_nmhc      ! normalised data
        CLOSE(21)
        OPEN(21,FILE=trim(emdir)//'bb_GFED_NH3_96_72.dat'//smonth,      &
     &    STATUS='OLD')
        READ(21,*) gfed_nh3       ! normalised data
        CLOSE(21)
        OPEN(21,FILE=trim(emdir)//'bb_GFED_SO2_96_72.dat'//smonth,      &
     &    STATUS='OLD')
        READ(21,*) gfed_so2       ! normalised data
        CLOSE(21)
      END IF

! bb_factor is a value for each year (1850-2100), relative to 1990,
! based on population in the tropics (1950-2050)/world(other years)
!      OPEN(21,FILE=trim(EMDIR)//'bb_pop_1850-2100.dat',STATUS='OLD')
!      READ(21,*) bb_factor
!      CLOSE(21)
!      bb_year=year-1850
!      bb_year=max(1,bb_year)
!      bb_year=min(251,bb_year)
!      burn=burn*bb_factor(bb_year)

      OPEN(22,FILE=trim(emdir)//'soilnox'//smonth//'_96_72.dat',        &
     &     STATUS='OLD')
      READ(22,*) soil
      CLOSE(22)
      OPEN(22,FILE=trim(emdir)//'mfch4wet'//smonth//'_96_72.dat',       &
     &     STATUS='OLD')
      READ(22,*) wetland
      CLOSE(22)
      OPEN(23,FILE=trim(emdir)//'nvoc_land'//smonth//'_96_72.dat',      &
     &     STATUS='OLD')
      READ(23,*) veg
      CLOSE(23)
      OPEN(23,FILE=trim(emdir)//'nvoc_ocean'//smonth//'_96_72.dat',     &
     &     STATUS='OLD')
      READ(23,*) ocean
      CLOSE(23)
      OPEN(23,FILE=trim(emdir)//'dms_ocean'//SMONTH//'_96_72.dat',      &
     &     STATUS='OLD')
      READ(23,*) OCEAN_DMS
      CLOSE(23)
!     OPEN(23,FILE=trim(emdir)//'isop'//smonth//'_96_72.dat',
!    &     STATUS='OLD')
!     READ(23,*) isop
!     CLOSE(23)
!     OPEN(23,FILE=trim(emdir)//'terp'//smonth//'_96_72.dat',
!    &     STATUS='OLD')
!     READ(23,*) terp
!     CLOSE(23)
! RADON data, kg(RN)/yr per grid square
      OPEN(24,FILE=trim(emdir)//'radon_96_72.dat',STATUS='OLD')
      READ(24,*) radon
      CLOSE(24)

! NOx       ! n.b. class(5) is non-seasonal shipping
!           ! n.b. burn is normalised, normalise others by ann. total
      IF (l_bu_selected) burn = gfed_nox
      emiss(i_no,:,:)=emiss(i_no,:,:)+                                  &
     &  (class(2,i_no)*burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+    &
     &  (class(3,i_no)/Annual_Land_NVOC)*                               &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_no)/Annual_Soil_NOx)*                                &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                 &
     &   na/(mnit*ysec)
! CO
      IF (l_bu_selected) burn = gfed_co
      emiss(i_co,:,:)=emiss(i_co,:,:)+                                  &
     &  (class(2,i_co)*burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+    &
     &  (class(3,i_co)/Annual_Land_NVOC)*                               &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_co)/Annual_Soil_NOx)*                                &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_co)/Annual_Ocean_NVOC)*                              &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mco*ysec)
! CH4  27 Tg/yr termites (=vegetation), 158 Tg/yr wetland
      IF (l_bu_selected) burn = gfed_ch4
      emiss(i_ch4,:,:)=emiss(i_ch4,:,:)+                                &
     &  (wetland(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)*wetl_fac+      &
     &   class(2,i_ch4)*                                                &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     & ((class(3,i_ch4)+27.0e9)/Annual_Land_NVOC)*                      &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_ch4)/Annual_Soil_NOx)*                               &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_ch4)/Annual_Ocean_NVOC)*                             &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mch4*ysec)
! SO2
      IF (l_bu_selected) burn = gfed_so2
      emiss(i_so2,:,:)=emiss(i_so2,:,:)+                                &
     &  (class(2,i_so2)*                                                &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_so2)/Annual_Land_NVOC)*                              &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_so2)/Annual_Soil_NOx)*                               &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_so2)/Annual_Ocean_NVOC)*                             &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(msul*ysec)
! NH3
      IF (l_bu_selected) burn = gfed_nh3
      emiss(i_nh3,:,:)=emiss(i_nh3,:,:)+                                &
     &  (class(2,i_nh3)*                                                &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_nh3)/Annual_Land_NVOC)*                              &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_nh3)/Annual_Soil_NOx)*                               &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_nh3)/Annual_DMS)*                                    &
     &   ocean_dms(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*            &
     &   na/(mnit*ysec)
!
      IF (l_bu_selected) burn = gfed_nmhc
! HCHO
      emiss(i_hcho,:,:)=emiss(i_hcho,:,:)+                              &
     &  (class(2,i_hcho)*                                               &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_hcho)/Annual_Land_NVOC)*                             &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_hcho)/Annual_Soil_NOx)*                              &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_hcho)/Annual_Ocean_NVOC)*                            &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mhcho*ysec)
! H2
      emiss(i_h2,:,:)=emiss(i_h2,:,:)+                                  &
     &  (class(2,i_h2)*                                                 &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_h2)/Annual_Land_NVOC)*                               &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_h2)/Annual_Soil_NOx)*                                &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_h2)/Annual_Ocean_NVOC)*                              &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mh2*ysec)
! C2H6
      emiss(i_c2h6,:,:)=emiss(i_c2h6,:,:)+                              &
     &  (class(2,i_c2h6)*                                               &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_c2h6)/Annual_Land_NVOC)*                             &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_c2h6)/Annual_Soil_NOx)*                              &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_c2h6)/Annual_Ocean_NVOC)*                            &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mc2h6*ysec)
! CH3CHO
      emiss(i_ch3cho,:,:)=emiss(i_ch3cho,:,:)+                          &
     &  (class(2,i_ch3cho)*                                             &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_ch3cho)/Annual_Land_NVOC)*                           &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_ch3cho)/Annual_Soil_NOx)*                            &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_ch3cho)/Annual_Ocean_NVOC)*                          &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mch3cho*ysec)
! NC4H10
      emiss(i_nc4h10,:,:)=emiss(i_nc4h10,:,:)+                          &
     &  (class(2,i_nc4h10)*                                             &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_nc4h10)/Annual_Land_NVOC)*                           &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_nc4h10)/Annual_Soil_NOx)*                            &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_nc4h10)/Annual_Ocean_NVOC)*                          &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mnc4h10*ysec)
! C2H4
      emiss(i_c2h4,:,:)=emiss(i_c2h4,:,:)+                              &
     &  (class(2,i_c2h4)*                                               &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_c2h4)/Annual_Land_NVOC)*                             &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_c2h4)/Annual_Soil_NOx)*                              &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_c2h4)/Annual_Ocean_NVOC)*                            &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mc2h4*ysec)
! C3H6
      emiss(i_c3h6,:,:)=emiss(i_c3h6,:,:)+                              &
     &  (class(2,i_c3h6)*                                               &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_c3h6)/Annual_Land_NVOC)*                             &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_c3h6)/Annual_Soil_NOx)*                              &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_c3h6)/Annual_Ocean_NVOC)*                            &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mc3h6*ysec)
! C3H8
      emiss(i_c3h8,:,:)=emiss(i_c3h8,:,:)+                              &
     &  (class(2,i_c3h8)*                                               &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_c3h8)/Annual_Land_NVOC)*                             &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_c3h8)/Annual_Soil_NOx)*                              &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_c3h8)/Annual_Ocean_NVOC)*                            &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mc3h8*ysec)
! Xylene
      emiss(i_oxyl,:,:)=emiss(i_oxyl,:,:)+                              &
     &  (class(2,i_oxyl)*                                               &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_oxyl)/Annual_Land_NVOC)*                             &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_oxyl)/Annual_Soil_NOx)*                              &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_oxyl)/Annual_Ocean_NVOC)*                            &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(moxyl*ysec)
! Toluene
      emiss(i_toluen,:,:)=emiss(i_toluen,:,:)+                          &
     &  (class(2,i_toluen)*                                             &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_toluen)/Annual_Land_NVOC)*                           &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_toluen)/Annual_Soil_NOx)*                            &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_toluen)/Annual_Ocean_NVOC)*                          &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mtoluen*ysec)
! DMS
      emiss(i_dms,:,:)=emiss(i_dms,:,:)+                                &
     &  (class(2,i_dms)*                                                &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_dms)/Annual_Land_NVOC)*                              &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_dms)/Annual_Soil_NOx)*                               &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_dms)/Annual_DMS)*                                    &
     &   ocean_dms(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)              &
     &   )*na/(msul*ysec)
! CH3OH
      emiss(i_ch3oh,:,:)=emiss(i_ch3oh,:,:)+                            &
     &  (class(2,i_ch3oh)*                                              &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_ch3oh)/Annual_Land_NVOC)*                            &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_ch3oh)/Annual_Soil_NOx)*                             &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_ch3oh)/Annual_Ocean_NVOC)*                           &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mch3oh*ysec)
! Acetone
      emiss(i_acetone,:,:)=emiss(i_acetone,:,:)+                        &
     &  (class(2,i_acetone)*                                            &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_acetone)/Annual_Land_NVOC)*                          &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_acetone)/Annual_Soil_NOx)*                           &
     &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(5,i_acetone)/Annual_Ocean_NVOC)*                         &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(macetone*ysec)
! C5H8
!     emiss(i_c5h8,:,:)=emiss(i_c5h8,:,:)+
!    &  (class(2,i_c5h8)*
!    &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+
!    &  (class(3,i_c5h8)/Annual_isop)*
!    &   isop(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+
!    &  (class(4,i_c5h8)/Annual_Soil_NOx)*
!    &   soil(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+
!    &  (class(5,i_c5h8)/Annual_Ocean_NVOC)*
!    &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*
!    &   na/(mc5h8*ysec)
! Radon
      emiss(i_rn222,:,:)=emiss(i_rn222,:,:)+                            &
     &  (class(2,i_rn222)*                                              &
     &   burn(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                  &
     &  (class(3,i_rn222)/Annual_Land_NVOC)*                            &
     &   veg(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                   &
     &  (class(4,i_rn222)/SUM(radon))*                                  &
     &   radon(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)+                 &
     &  (class(5,i_rn222)/Annual_Ocean_NVOC)*                           &
     &   ocean(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1))*                &
     &   na/(mrn222*ysec)

      END SUBROUTINE EMUPDT
#endif
