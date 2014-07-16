#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EMCALC_PI(emiss,month,year)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : CALCULATE EMISSIONS
!-
!-   Inputs  : MONTH
!-   Outputs : EMISS
!-   Controls:
!-
!
! Current Code Owner: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.5    23/12/93  Created.  W.J. Collins
!  6.1    14/10/04  Revised and improved emissions. C.E. Johnson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!
!-
!VVV  v4.5  EMCALC_PI 11/12/00 - CLASS in kg, not g.
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: month
      INTEGER, INTENT(IN) :: year
      REAL, DIMENSION(nc,nlnpe,nlpe), INTENT(OUT) :: emiss

      INTEGER :: k
      INTEGER :: ic
      INTEGER :: iyr1
      INTEGER :: iyr2
      REAL, DIMENSION(nlnpe,nlpe) :: isopre
      REAL, DIMENSION(5,nc) :: class

      class=0.
! Note that termites and wetland CH4 emissions are specified in
! EMUPDT and EMREAD respectively.
! Paddy and animals are assumed to be part of the anthropogenic.

!
!     -Pre-Industrial----EMISSION CLASSES-(Tg/yr)----------------------
!
! PI  1850 (from emiss_pi.pro)
!                            Anth  ,Biomass,   Veg ,  Soil , Oceans
      class(:,     i_no)= (/   0.00,   1.02,   0.00,   3.50,   6.90 /)
      class(:,     i_co)= (/   0.00,  50.70, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/   0.00,   2.36,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/   0.00,   2.00,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.00,   0.06,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   0.00,   0.33,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   0.00,   0.29,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/   0.00,   0.16,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   0.00,   0.69,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   0.00,   0.33,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   0.00,   0.09,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   0.00,   0.38, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   0.00,   0.04,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/   0.00,   0.44,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   0.00,   0.30,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/   0.00,   0.59,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/   0.00,   0.14,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   1.00,  24.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)

      class=class*1e9
      WRITE(6,*) '**** PRE-INDUSTRIAL **** ',YEAR,' ****'
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
        IF(k/=i_o3.AND.k/=i_hno3) emiss(k,:,:)=emiss(k,:,:)/lmolec
      END DO

      END SUBROUTINE EMCALC_PI
#endif
