#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ To calculate emissions for STOCHEM from month and year
!
! Subroutine interface
      SUBROUTINE EMCALC_MFR(emiss,month,year)
!
! Description
!  Emissions for Dentener/Stevenson scenarios for 2000, 2030
!
! Method
!  developed using wave routine: emiss_mfr.pro
!
! Current Code Owner: Colin Johnson
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  6.2    01/03/06  Created for ACCENT IPCC simulations. W.J. Collins
!
! Code Description:
!   Language: FORTRAN 90 + common extensions.
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!----------------------------------------------------------------------
! Subroutine arguments
      INTEGER, INTENT(IN) :: month
      INTEGER, INTENT(IN) :: year
      REAL, DIMENSION(nc,nlnpe,nlpe), INTENT(OUT) :: emiss

! Local
      INTEGER :: k
      INTEGER :: ic
      INTEGER :: iyr1
      INTEGER :: iyr2
      INTEGER :: yearx
      REAL, DIMENSION(5,nc,2)       :: yclass
      REAL, DIMENSION(5,nc)         :: class
      CHARACTER(72) :: cmessage

      class = 0.0

! Note that termites and wetland CH4 emissions are specified in
! EMUPDT and EMREAD respectively.
! Paddy and animals are assumed to be part of the anthropogenic.

!     -------------------EMISSION classes-(Tg/yr)----------------------
! MFR     2000
!                            Anth  ,Biomass,  Veg  , Soil  , Oceans
      class(:,     i_no)= (/  27.80,  10.20,   0.00,   5.60,   0.00 /)
      class(:,     i_co)= (/ 470.00, 507.00, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 300.50,  23.60,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  20.00,  20.00,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.77,   0.57,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   7.20,   3.32,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.45,   2.89,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  58.97,   1.56,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   7.97,   6.92,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   3.68,   3.32,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   8.42,   0.95,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   4.90,   3.84, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   2.91,   0.38,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/  10.26,   4.41,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   8.27,   3.03,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  54.20,   5.90,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  54.20,   1.40,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   1.00,  24.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,1) = class(:,:) * 1.0e9

! MFR     2030
!                            Anth   , Biomass ,   Veg   ,  Soil   , Ocea
      class(:,     i_no)= (/  13.10,  10.20,   0.00,   5.60,   0.00 /)
      class(:,     i_co)= (/ 221.70, 507.00, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 290.00,  23.60,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/   9.42,  20.00,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.48,   0.57,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   4.55,   3.32,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   1.55,   2.89,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  37.28,   1.56,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   5.03,   6.92,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   2.32,   3.32,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   5.33,   0.95,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   3.10,   3.84, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   1.84,   0.38,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/   6.49,   4.41,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   5.23,   3.03,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  65.00,   4.86,   0.00,   2.40,   8.30 /)
      class(:,    i_so2)= (/  16.50,   1.40,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   1.00,  24.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,2) = class(:,:) * 1.0e9

! Set interpolation year
      IF (l_emiss_current) THEN
        yearx = year
      ELSE
        yearx = emiss_year
        WRITE(6,*) 'L_Emiss_Year is False, using emissions from: ',     &
     &    emiss_year
      END IF

! Find right year
      IF (yearx == 2000) THEN
        class(:,:)=yclass(:,:,1)
      ELSE IF (yearx == 2030) THEN
        class(:,:)=yclass(:,:,2)
      ELSE
        WRITE(6,*) 'Wrong year selected for MFR emissions, yearx: ',year
        cmessage='Wrong year selected for MFR emissions'
! DEPENDS ON: ereport
        CALL EREPORT('Emcalc_MFR',1,cmessage)
      ENDIF

      WRITE(6,*) '**** MFR **** ',yearx,' ****'

! Get EMISS in units (molecules s^-1 per grid square)
! DEPENDS ON: emread
      CALL EMREAD(emiss,year,class)

! Add seasonal variations in biomass burning, soil
! & vegetation emissions
! DEPENDS ON: emupdt
      CALL EMUPDT(emiss,class,month,year)

! Convert emissions to (molecules s^-1) per grid square /
! (molecules/cell)
      DO k=1,nc                        !
! Don't scale O3 and HNO3 (will do in STRATCALC)
        IF (k /= i_o3 .AND. k /= i_hno3)                                &
     &    emiss(k,:,:) = emiss(k,:,:) / lmolec
      END DO

      END SUBROUTINE EMCALC_MFR
#endif
