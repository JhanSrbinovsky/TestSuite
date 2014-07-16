#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EMCALC_B1(emiss,month,year)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Calculate SRES B1 scenario emissions
!-
!-   Inputs  : MONTH
!-   Outputs : EMISS
!-   Controls:
!
! Current Code Owner: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.3    23/12/93  Created.  W.J. Collins
!  4.5    04/10/99  Added year dependence: YCLASS arrays. C.E. Johnson
!  5.5    13/01/04  Calculates correct emissions when l_emiss_current
!                   is false. C.E. Johnson
!  6.2    01/03/06  Removed calculation of isoprene emission as now done
!                   interactively.  M.G. Sanderson
!
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
      INTEGER :: yearx
      REAL, DIMENSION(5,nc,198:211) :: yclass
      REAL, DIMENSION(5,nc)         :: class
      CHARACTER*72 :: cmessage

      class = 0.0
! Note that termites and wetland CH4 emissions are specified in
! EMUPDT and EMREAD respectively.
! Paddy and animals are assumed to be part of the anthropogenic.
!
! SRES B1    1990
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  31.00,   7.10,   0.00,   5.60,   6.68 /)
      class(:,     i_co)= (/ 650.00, 500.00, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 310.00,  40.00,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  20.00,  20.00,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.92,   0.60,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   8.64,   3.51,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.94,   3.06,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  70.79,   1.66,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   9.56,   7.32,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   4.41,   3.51,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/  10.11,   1.00,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   5.88,   4.06, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   3.49,   0.40,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/  12.32,   4.66,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   9.93,   3.21,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  39.40,   3.50,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  70.90,   2.20,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,199) = class(:,:) * 1.0e9
!
! SRES B1    2000
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  32.00,   8.46,   0.00,   5.60,   6.90 /)
      class(:,     i_co)= (/ 648.00, 596.09, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 323.00,  47.69,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  20.65,  23.84,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.93,   0.72,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   8.77,   4.19,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.98,   3.65,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  71.81,   1.97,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   9.70,   8.73,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   4.48,   4.19,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/  10.26,   1.20,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   5.97,   4.84, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   3.54,   0.48,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/  12.50,   5.56,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  10.07,   3.83,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  40.67,   4.17,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  69.00,   2.62,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,200) = class(:,:) * 1.0e9
!
! SRES B1    2010
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  36.00,   9.67,   0.00,   5.60,   7.76 /)
      class(:,     i_co)= (/ 560.00, 681.25, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 349.00,  54.50,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  23.23,  27.25,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.93,   0.82,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   8.77,   4.78,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.98,   4.17,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  71.81,   2.25,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   9.70,   9.98,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   4.48,   4.78,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/  10.26,   1.37,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   5.97,   5.53, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   3.54,   0.55,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/  12.50,   6.35,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  10.07,   4.37,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  45.75,   4.77,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  73.90,   3.00,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,201) = class(:,:) * 1.0e9
!
! SRES B1    2020
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  40.00,  10.81,   0.00,   5.60,   8.62 /)
      class(:,     i_co)= (/ 522.00, 760.94, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 377.00,  60.88,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  25.81,  30.44,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.93,   0.92,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   8.70,   5.34,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.96,   4.66,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  71.30,   2.52,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   9.63,  11.14,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   4.44,   5.34,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/  10.19,   1.53,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   5.93,   6.18, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   3.52,   0.61,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/  12.41,   7.10,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  10.00,   4.88,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  50.84,   5.33,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  74.60,   3.35,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,202) = class(:,:) * 1.0e9
!
! SRES B1    2030
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  42.00,  11.69,   0.00,   5.60,   9.06 /)
      class(:,     i_co)= (/ 374.00, 823.44, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 385.00,  65.88,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  27.10,  32.94,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.87,   0.99,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   8.14,   5.78,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.77,   5.04,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  66.71,   2.73,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   9.01,  12.06,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   4.16,   5.78,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   9.53,   1.65,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   5.54,   6.69, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   3.29,   0.66,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/  11.61,   7.68,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   9.36,   5.29,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  53.38,   5.76,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  78.20,   3.62,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,203) = class(:,:) * 1.0e9
!
! SRES B1    2040
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  43.00,  12.23,   0.00,   5.60,   9.27 /)
      class(:,     i_co)= (/ 302.00, 860.94, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 381.00,  68.88,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  27.74,  34.44,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.81,   1.04,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   7.65,   6.04,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.60,   5.27,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  62.64,   2.85,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   8.46,  12.61,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   3.90,   6.04,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   8.95,   1.73,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   5.21,   6.99, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   3.09,   0.69,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/  10.90,   8.03,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   8.79,   5.53,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  54.65,   6.03,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  78.50,   3.79,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,204) = class(:,:) * 1.0e9
!
! SRES B1    2050
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  39.00,  12.46,   0.00,   5.60,   8.41 /)
      class(:,     i_co)= (/ 242.00, 877.34, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 359.00,  70.19,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  25.16,  35.09,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.77,   1.06,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   7.21,   6.16,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.46,   5.37,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  59.07,   2.90,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   7.98,  12.85,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   3.68,   6.16,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   8.44,   1.76,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   4.91,   7.13, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   2.92,   0.70,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/  10.28,   8.18,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   8.29,   5.63,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  49.57,   6.14,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  68.90,   3.86,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,205) = class(:,:) * 1.0e9
!
! SRES B1    2060
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  34.00,  12.37,   0.00,   5.60,   7.33 /)
      class(:,     i_co)= (/ 230.00, 871.09, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 342.00,  69.69,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  21.94,  34.84,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.73,   1.05,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   6.90,   6.12,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.35,   5.33,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  56.53,   2.88,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   7.63,  12.76,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   3.52,   6.12,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   8.08,   1.75,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   4.70,   7.08, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   2.79,   0.70,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/   9.84,   8.13,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   7.93,   5.59,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  43.21,   6.10,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  55.80,   3.83,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,206) = class(:,:) * 1.0e9
!
! SRES B1    2070
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  30.00,  11.86,   0.00,   5.60,   6.47 /)
      class(:,     i_co)= (/ 227.00, 835.16, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 324.00,  66.81,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  19.35,  33.41,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.68,   1.01,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   6.40,   5.86,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.18,   5.11,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  52.45,   2.76,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   7.08,  12.23,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   3.27,   5.86,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   7.49,   1.68,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   4.36,   6.79, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   2.59,   0.67,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/   9.13,   7.79,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   7.36,   5.36,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  38.13,   5.85,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  44.30,   3.67,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,207) = class(:,:) * 1.0e9
!
! SRES B1    2080
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  26.00,  11.49,   0.00,   5.60,   5.61 /)
      class(:,     i_co)= (/ 197.00, 809.38, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 293.00,  64.75,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  16.77,  32.38,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.65,   0.97,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   6.15,   5.68,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.10,   4.95,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  50.42,   2.68,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   6.81,  11.85,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   3.14,   5.68,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   7.20,   1.62,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   4.19,   6.58, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   2.49,   0.65,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/   8.77,   7.55,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   7.07,   5.20,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  33.05,   5.67,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  36.10,   3.56,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,208) = class(:,:) * 1.0e9
!
! SRES B1    2090
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  22.00,  10.72,   0.00,   5.60,   4.74 /)
      class(:,     i_co)= (/ 170.00, 754.69, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 266.00,  60.38,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  14.19,  30.19,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.63,   0.91,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   5.97,   5.30,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   2.03,   4.62,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  48.89,   2.50,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   6.60,  11.05,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   3.05,   5.30,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   6.98,   1.51,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   4.06,   6.13, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   2.41,   0.61,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/   8.51,   7.04,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   6.86,   4.84,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  27.96,   5.28,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  29.80,   3.32,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,209) = class(:,:) * 1.0e9
!
! SRES B1    2100
!                            Anth  , Biomass ,  Veg ,  Soil, Oceans
      class(:,     i_no)= (/  19.00,   9.84,   0.00,   5.60,   4.10 /)
      class(:,     i_co)= (/ 134.00, 692.97, 130.00,   0.00,  40.00 /)
      class(:,    i_ch4)= (/ 236.00,  55.44,   0.00,   0.00,  10.00 /)
      class(:,     i_h2)= (/  12.26,  27.72,   0.00,   4.00,   4.00 /)
      class(:,   i_hcho)= (/   0.58,   0.83,   0.00,   0.00,   0.00 /)
      class(:,   i_c2h6)= (/   5.41,   4.87,   1.60,   0.30,   1.50 /)
      class(:, i_ch3cho)= (/   1.84,   4.24,   0.00,   0.00,   0.00 /)
      class(:, i_nc4h10)= (/  44.31,   2.29,  63.10,   0.10,   8.30 /)
      class(:,   i_c2h4)= (/   5.98,  10.15,  11.80,   3.20,   8.30 /)
      class(:,   i_c3h6)= (/   2.76,   4.87,   7.70,   0.20,   1.40 /)
      class(:,   i_c3h8)= (/   6.33,   1.39,   1.60,   0.90,   6.40 /)
      class(:,  i_ch3oh)= (/   3.68,   5.63, 100.00,   0.00,   0.00 /)
      class(:,i_acetone)= (/   2.19,   0.56,  40.00,   0.00,   0.00 /)
      class(:, i_toluen)= (/   7.71,   6.46,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   6.21,   4.45,   0.00,   0.00,   0.00 /)
      class(:,   i_c5h8)= (/   0.00,   0.00, 450.00,   0.00,   0.00 /)
      class(:,    i_nh3)= (/  24.15,   4.85,   0.00,   2.50,   8.20 /)
      class(:,    i_so2)= (/  24.90,   3.05,   0.00,   0.00,   0.00 /)
      class(:,    i_dms)= (/   0.00,   0.00,   0.00,   0.90,  15.00 /)
      class(:,  i_rn222)= (/   0.00,   0.00,   0.00,15.0e-9,   0.00 /)
      yclass(:,:,210) = class(:,:) * 1.0e9

      yclass(:,:,211) = yclass(:,:,210)
      yclass(:,:,198) = yclass(:,:,199)

! Set interpolation year
      IF (l_emiss_current) THEN
        yearx = year
      ELSE
        yearx = emiss_year
        WRITE(6,*) 'L_Emiss_Year is False, using emissions from: ',     &
     &    emiss_year
      END IF

! Interpolate to find class
      iyr1 = yearx / 10
      IF (l_emiss_current) THEN
        iyr2 = iyr1
      ELSE
        iyr2 = iyr1 + 1
      END IF
      IF (iyr2 > 211) THEN
        class(:,:) = yclass(:,:,210)
        cmessage = 'STOCHEM emissions have wrong date, using 2100'
        WRITE(6,*) cmessage
        WRITE(6,*) iyr1,iyr2,yearx
! DEPENDS ON: ereport
        CALL EREPORT('Emcalc_B1',-1,cmessage)
        WRITE(6,*) '**** SRES B1 **** ',2100,' ****'
      ELSE IF (iyr1 < 198) THEN
        class(:,:) = yclass(:,:,198)
        cmessage = 'STOCHEM emissions have wrong date, using 1990'
        WRITE(6,*) cmessage
        WRITE(6,*) iyr1,iyr2,yearx
! DEPENDS ON: ereport
        CALL EREPORT('EMCALC_B1',-1,cmessage)
        WRITE(6,*) '**** SRES B1 **** ',1990,' ****'
      ELSE
        class(:,:) = yclass(:,:,iyr1) + REAL(yearx-iyr1*10) *           &
     &    (yclass(:,:,iyr2)-yclass(:,:,iyr1)) / 10.0
        WRITE(6,*) '**** SRES B1 **** ',YEAR,' ****'
      END IF

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

      END SUBROUTINE EMCALC_B1
#endif
