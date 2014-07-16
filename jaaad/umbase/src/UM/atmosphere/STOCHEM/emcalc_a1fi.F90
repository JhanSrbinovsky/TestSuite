#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EMCALC_A1FI(emiss,month,year)
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
!  3.3    23/12/93  Created.  W.J. Collins
!  4.5    04/10/99  Added year dependence: YCLASS arrays. C.E. Johnson
!  5.5    13/01/04  Calculates correct emissions when l_emiss_current
!                   is false. C.E. Johnson
!  6.2    01/03/06  Added missing emissions. Removed isoprene emission
!                   as now calculated interactively.  M.G. Sanderson
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
      CHARACTER(72) :: cmessage

      class = 0.0
! Note that termites and wetland CH4 emissions are specified in
! EMUPDT and EMREAD respectively.
! Paddy and animals are assumed to be part of the anthropogenic.

!     -1990-A1F1---------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/23.6  ,     7.4,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/386.8 ,   492.2, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/259.4 ,    50.6,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/1.4   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/6.0   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/3.0   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/54.2  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/6.9   ,     8.5,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/9.4   ,     9.6,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/6.2   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/7.2   ,     8.9, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/2.5   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,     0.0, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_so2)=    (/70.9  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  12.32,   4.66,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/   9.93,   3.21,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,199)=class(:,:)*1e9

!     -2000-A1F1---------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/24.4  ,     7.6,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/386.0 ,   491.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/270.3 ,    52.7,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/1.4   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/6.1   ,     4.5,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/3.1   ,     4.9,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/54.9  ,     2.7,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/7.0   ,     8.6,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/9.5   ,     9.7,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/6.3   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/7.3   ,     9.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/2.6   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/69.0  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  12.50,   5.58,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  10.07,   3.84,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,200)=class(:,:)*1e9

!
!     -2010-A1F1---------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/30.5  ,     9.5,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/448.9 ,   571.1, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/300.4 ,    58.6,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/1.7   ,     1.8,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/7.2   ,     5.3,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/3.6   ,     5.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/64.7  ,     3.1,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/8.2   ,    10.1,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/11.2  ,    11.4,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/7.5   ,     1.3,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/8.5   ,    10.6, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/3.0   ,     0.8,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/80.8  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  14.71,   6.63,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  11.86,   4.56,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,201)=class(:,:)*1e9

!
!     -2020-A1F1---------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/43.6  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/734.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/408.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.1   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/9.1   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/4.5   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/81.9  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/10.4  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/14.2  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/9.5   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/10.8  ,     8.8, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/3.8   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/86.9  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  17.02,   7.68,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  13.71,   5.29,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,202)=class(:,:)*1e9

!
!     -2030-A2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/53.7  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/913.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/466.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.4   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/10.5  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/5.2   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/94.4  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h6)=   (/12.0  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/16.3  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/10.9  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/12.4  ,     8.8, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/4.3   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/96.1  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  18.97,   8.73,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  15.29,   6.01,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,203)=class(:,:)*1e9

!
!     -2040-A2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/58.7  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/994.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/520.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.6   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/11.2  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/5.6   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/100.6 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/12.8  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/17.4  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/11.6  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/13.2  ,     8.8, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/4.6   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/94.0  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  22.69,  10.24,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  18.29,   7.05,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,204)=class(:,:)*1e9

!
!     -2050-A2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/64.4  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1086.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/581.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.7   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/12.0  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/6.0   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/108.0 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/13.7  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/18.7  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/12.5  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/14.2  ,     8.8, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/5.0   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/80.5  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  28.54,  10.86,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  23.00,   7.47,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,205)=class(:,:)*1e9

!
!     -2060-A2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/68.5  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1196.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/634.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.9   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/12.7  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/6.3   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/114.3 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/14.5  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/19.7  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/13.2  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/15.0  ,     8.8, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/5.3   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/56.3  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  31.99,  11.70,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  25.79,   8.05,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,206)=class(:,:)*1e9

!
!     -2070-A2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/73.0  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1316.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/692.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.1   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/13.5  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/6.8   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/122.0 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/15.5  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/21.1  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/14.1  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/16.0  ,     8.8, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/5.6   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/42.6  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  35.89,  12.23,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  28.93,   8.42,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,207)=class(:,:)*1e9

!
!     -2080-A2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/80.2  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1485.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/750.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.4   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/14.9  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/7.5   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/134.7 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/17.1  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/23.3  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/15.6  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/17.7  ,     8.8, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/6.2   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/39.4  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  39.79,  13.48,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  32.07,   9.28,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,208)=class(:,:)*1e9

!
!     -2090-A2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/90.8  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1718.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/808.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.9   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/17.0  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/8.5   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/152.9 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/19.4  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/26.4  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/17.7  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/20.1  ,     8.8, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/7.0   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/39.8  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  38.55,  14.38,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  31.07,   9.89,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,209)=class(:,:)*1e9

!
!     -2100-A2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/102.6 ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1984.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/873.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/4.4   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/19.3  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/9.6   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/173.8 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/22.0  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/30.0  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/20.1  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/22.9  ,     8.8, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/8.0   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/40.1  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/  0.0 ,  0.00,    0.00,  0.90, 15.00 /)
      class(:, i_toluen)= (/  37.22,  14.68,   0.00,   0.00,   0.00 /)
      class(:,   i_oxyl)= (/  30.00,  10.10,   0.00,   0.00,   0.00 /)
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,210)=class(:,:)*1e9

      yclass(:,:,211) = yclass(:,:,210)
      yclass(:,:,198) = yclass(:,:,199)

! Set interpolation year
      IF (L_emiss_current) THEN
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
      if (iyr2 > 211) THEN
        class(:,:)=yclass(:,:,210)
        cmessage = 'STOCHEM emissions have wrong date, using 2100'
        WRITE(6,*) cmessage
        WRITE(6,*) iyr1,iyr2,yearx
! DEPENDS ON: ereport
        CALL EREPORT('Emcalc_A2',-1,cmessage)
        WRITE(6,*) '**** SRES A2 **** ',2100,' ****'
      ELSE IF (iyr1 < 198) THEN
        class(:,:)=yclass(:,:,198)
        cmessage = 'STOCHEM emissions have wrong date, using 1990'
        WRITE(6,*) cmessage
        WRITE(6,*) iyr1,iyr2,yearx
! DEPENDS ON: ereport
        CALL EREPORT('Emcalc_A2',-1,cmessage)
        WRITE(6,*) '**** SRES A2 **** ',1990,' ****'
      ELSE
        class(:,:)=yclass(:,:,iyr1)+REAL(yearx-iyr1*10)*                &
     &    (yclass(:,:,iyr2)-yclass(:,:,iyr1))/10.0
      END IF

      WRITE(6,*) '**** SRES A1FI **** ',yearx
! Get EMISS in units (molecules s^-1 per grid square)
! DEPENDS ON: emread
      CALL EMREAD(emiss,year,class)
! Add seasonal variations in biomass burning, soil
! & vegetation emissions
! DEPENDS ON: emupdt
      CALL EMUPDT(emiss,class,month,year)
! Convert emissions to (molecules s^-1) per grid square /
! (molecules/cell)

      DO k=1,nc
! Don't scale O3 and HNO3 (will do in STRATCALC)
        IF(k/=i_o3.AND.k/=i_hno3) emiss(k,:,:)=emiss(k,:,:)/lmolec
      END DO

      END SUBROUTINE EMCALC_A1FI
#endif
