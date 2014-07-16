#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EMCALC_A2(emiss,month,year)
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
!  6.2    10/03/06  Added emissions of aromatics and DMS. M.G. Sanderson
!
!VVV  v4.5  EMCALC_A2 11/12/00 - class in kg, not g.
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
      REAL, DIMENSION(nlnpe,nlpe)   :: isopre
      REAL, DIMENSION(5,nc,198:211) :: yclass
      REAL, DIMENSION(5,nc)         :: class
      CHARACTER(len=72) :: cmessage

      class = 0.0

! Note that termites and wetland CH4 emissions are specified in
! EMUPDT and EMREAD respectively.
! Paddy and animals are assumed to be part of the anthropogenic.

!     -1990-A2-----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/23.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/379.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/255.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/1.4   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/6.1   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/3.1   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/54.2  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/6.9   ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/9.5   ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/6.3   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/7.2   ,    13.0, 100.0,  0.0,  0.0 /) ! CH3OH
      class(:,i_acetone)=(/2.5   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/70.9  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/8.93  ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/9.01  ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.e-9, 0.0 /) ! RN222
      yclass(:,:,199)=class(:,:)*1.0e9
!
!     -2000-A2-----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/24.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/377.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/268.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/1.5   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/6.2   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/3.2   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/55.0  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/7.2   ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/9.7  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/6.4   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/7.5   ,    13.0, 100.0,  0.0,  0.0 /) ! CH3OH
      class(:,i_acetone)=(/2.6   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/69.0  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/8.93  ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/9.01  ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,200)=class(:,:)*1.0e9
!
!     -2010-A2-----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/31.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/477.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/315.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/1.8   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/7.3   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/4.0   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/60.7  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/8.7   ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/11.6  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/7.1   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/9.1   ,    13.0, 100.0,  0.0,  0.0 /) ! CH3OH
      class(:,i_acetone)=(/2.9   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/74.7  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/9.82  ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/9.90  ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,201)=class(:,:)*1.0e9
!
!     -2020-A2-----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/42.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/575.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/369.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.3   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/9.1   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/5.3   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/70.5  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/11.3  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/14.9  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/8.4   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/11.9  ,    13.0, 100.0,  0.0,  0.0 /) ! CH3OH
      class(:,i_acetone)=(/3.5   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/99.5  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/11.34 ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/11.44 ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,202)=class(:,:)*1.0e9
!
!     -2030-A2-----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/53.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/759.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/431.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.7   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/10.8  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/6.6   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/79.9  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/13.9  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/18.1  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/9.6  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/14.5  ,    13.0, 100.0,  0.0,  0.0 /) ! CH3OH
      class(:,i_acetone)=(/4.0   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/112.5 ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/12.79 ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/12.91 ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,203)=class(:,:)*1.0e9
!
!     -2040-A2-----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/58.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/844.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/487.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.0   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/11.7  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/7.3   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/84.8 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/15.2  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/19.7  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/10.2  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/15.9  ,    13.0, 100.0,  0.0,  0.0 /) ! CH3OH
      class(:,i_acetone)=(/4.3   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/109.0 ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/13.55 ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/13.67 ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,204)=class(:,:)*1.0e9
!
!     -2050-A2-----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/63.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/928.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/543.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.2   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/12.6  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/7.9   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/89.3 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/16.4  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/21.2  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/10.8  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/17.1  ,    13.0, 100.0,  0.0,  0.0 /) ! CH3OH
      class(:,i_acetone)=(/4.5   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/105.4 ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/14.25 ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/14.38 ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,205)=class(:,:)*1.0e9
!
!     -2060-A2-----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/67.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1045.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/599.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.5   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/13.6  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/8.7   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/94.6 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/17.8  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/23.0  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/11.5  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/18.7  ,    13.0, 100.0,  0.0,  0.0 /) ! CH3OH
      class(:,i_acetone)=(/4.8   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/89.6  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/15.07 ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/15.21 ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,206)=class(:,:)*1.0e9
!
!     -2070-A2-----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/72.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1162.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/656.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.8   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/14.5  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/9.4   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/99.9 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/19.3  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/24.7  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/12.2  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/20.2  ,    13.0, 100.0,  0.0,  0.0 /) ! CH3OH
      class(:,i_acetone)=(/5.1   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/73.7  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/15.9  ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/16.04 ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN222
      yclass(:,:,207)=class(:,:)*1.0e9
!
!     -2080-A2-----------EMISSION classES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/79.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1342.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/715.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/4.3   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/16.3  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/10.7  ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/109.7 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/21.9  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/28.0  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/13.5  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/22.9  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/5.7   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/64.7  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/17.42 ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/17.57 ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN222
      yclass(:,:,208)=class(:,:)*1.0e9
!
!     -2090-A2-----------EMISSION classES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/90.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1584.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/774.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/5.0   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/18.9  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/12.7  ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/123.6 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/25.7  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/32.7  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/15.3  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/26.8  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/6.5   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/62.5  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/19.57 ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/19.74 ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN222
      yclass(:,:,209)=class(:,:)*1.0e9

!
!     -2100-A2-----------EMISSION classES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/101.9 ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1826.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/834.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/5.7   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/21.4  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/14.5  ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/137.1 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/29.3  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/37.1  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/17.0  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/30.7  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/7.2   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/60.3  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_toluen)= (/21.66 ,    2.27,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/21.85 ,    0.46,   0.0,  0.0,  0.0 /) ! OXYL
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN222
      yclass(:,:,210)=class(:,:)*1.0e9

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
        CALL EREPORT('Emcalc_A2',-1,cmessage)
        WRITE(6,*) '**** SRES A2 **** ',2100,' ****'
      ELSE IF (iyr1 < 198) THEN
        class(:,:) = yclass(:,:,198)
        cmessage = 'STOCHEM emissions have wrong date, using 1990'
        WRITE(6,*) cmessage
        WRITE(6,*) iyr1,iyr2,yearx
! DEPENDS ON: ereport
        CALL EREPORT('Emcalc_A2',-1,cmessage)
        WRITE(6,*) '**** SRES A2 **** ',1990,' ****'
      ELSE
        class(:,:) = yclass(:,:,iyr1) + REAL(yearx-iyr1*10) *           &
     &    (yclass(:,:,iyr2) - yclass(:,:,iyr1)) / 10.0
        WRITE(6,*) '**** SRES A2 **** ',yearx,' ****'
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

      END SUBROUTINE EMCALC_A2
#endif
