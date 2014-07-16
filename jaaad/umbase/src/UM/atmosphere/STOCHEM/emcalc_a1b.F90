#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EMCALC_A1B(emiss,month,year)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : CALCULATE EMISSIONS SRES A1B (AIM)
!-
!-   Inputs  : MONTH
!-   Outputs : EMISS
!-   Controls:
!
! Current Code Owner: W.J. Collins
!
! History:
! Version   Date                    Comment
!  6.1  20/08/04  New UM deck.  C.E. Johnson
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
      REAL, DIMENSION(nlnpe,nlpe)   :: isopre
      REAL, DIMENSION(5,nc,198:211) :: yclass
      REAL, DIMENSION(5,nc)         :: class
      CHARACTER (len=72) :: cmessage

      class = 0.0
! Note that termites and wetland CH4 emissions are specified in
! EMUPDT and EMREAD respectively.
! Paddy and animals are assumed to be part of the anthropogenic.

!     -1990-A1B----------EMISSION classes-(Tg/yr)----------------------
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
      CLASS(:,I_CH3OH)=  (/7.2   ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/2.5   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/70.9  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,199)=class(:,:)*1.0e9

!
!     -2000-A1B----------EMISSION classes-(Tg/yr)----------------------
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
      CLASS(:,I_CH3OH)=  (/7.5   ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/2.6   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/69.0  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,200)=class(:,:)*1.0e9

!
!     -2010-A1B----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/31.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/502.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/318.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.2   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/9.0   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/5.3   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/70.1  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/11.2  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/14.8  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/8.3   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/11.7  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/3.4   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/74.7  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,201)=class(:,:)*1.0e9

!
!     -2020-A1B----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/38.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/609.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/411.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/4.1   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/15.7  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/10.3  ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/106.1 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/20.9  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/26.8  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/13.0  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/21.9  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/5.5   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/99.5  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,202)=class(:,:)*1.0e9

!
!     -2030-A1B----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/42.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/609.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/411.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/4.1   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/15.7  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/10.2  ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/106.1 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/20.9  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/26.8  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/13.0  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/21.9  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/5.5   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/105.4 ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,203)=class(:,:)*1.0e9

!
!     -2040-A1B----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/41.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/660.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/403.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/4.2   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/16.1  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/10.6  ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/108.5 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/21.6  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/27.6  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/13.3 ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/22.6  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/5.6   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/112.5 ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,204)=class(:,:)*1.0e9

!
!     -2050-A1B----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/40.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/714.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/397.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/4.4   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/16.7  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/11.0  ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/111.4 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/22.4  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/28.6  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/13.7  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/23.4  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/5.8   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/109.0 ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,205)=class(:,:)*1.0e9

!
!     -2060-A1B----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/38.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/745.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/355.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/4.5   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/17.0  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/11.3  ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/113.4 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/22.9  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/29.2  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/13.9  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/24.0  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/5.9   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/105.4 ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,206)=class(:,:)*1.0e9

!
!     -2070-A1B----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/36.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/776.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/318.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/4.6   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/17.4  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/11.6  ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/115.5 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/23.5  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/29.9  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/14.2  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/24.5  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/6.0   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/89.6  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,207)=class(:,:)*1.0e9

!
!     -2080-A1B----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/35.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/857.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/286.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/4.2   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/15.9  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/10.4  ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/107.3 ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/21.3  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/27.2  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/13.1  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/22.2  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/5.5   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/73.7  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,208)=class(:,:)*1.0e9

!
!     -2090-A1B----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/33.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/999.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/259.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.3   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/12.8  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/8.1   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/90.5  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/16.7  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/21.6  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/11.0  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/17.5  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/4.6   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/64.7  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,209)=class(:,:)*1.0e9

!
!     -2100-A1B----------EMISSION classes-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/32.9  ,     7.1,   0.0,  5.6,  0.0 /) ! NO
      class(:,i_co)=     (/1163.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/234.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.6   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/10.2  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/6.1   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/76.2  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/12.9  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/16.8  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/9.1   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      CLASS(:,I_CH3OH)=  (/13.5  ,    13.0, 100.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/3.8   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/62.5  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,210) = class(:,:)*1.0e9

      yclass(:,:,211) = yclass(:,:,210)
      yclass(:,:,198) = yclass(:,:,199)

! Mod to add 425 Tg to anthropogenic emissions for specified year,month
!      IF (YEAR == 1990 .AND. MONTH == 1) THEN
!        yclass(1,I_CO,199) = yclass(1,I_CO,199) + 425.0*1.0E12
!        WRITE(6,*) 'YEAR/MONTH ',year,'/',month,' 500 Tg ',
!     &    'added to anthropogenic EXTRA'
!      ELSE   ! to put back to original
!        yclass(1,I_CO,199) = 393.0*1.0E12
!      END IF

! Interpolate to find class
      IF (.NOT. l_emiss_current) THEN
! Use emissions from a specified year
        WRITE(6,*) ' ***** Current emissions NOT being used in STOCHEM'
        WRITE(6,*) ' ***** Using emissions from year: ',Emiss_Year
        iyr1 = emiss_year / 10
        yearx = emiss_year
      ELSE
        iyr1 = year / 10
        yearx = year
      ENDIF
      iyr2 = iyr1 + 1
      IF (iyr2 > 211) THEN
        class(:,:) = yclass(:,:,210)
        cmessage = 'STOCHEM emissions have wrong date, using 2100'
        WRITE(6,*) cmessage
        WRITE(6,*) iyr1,iyr2,year
! DEPENDS ON: ereport
        CALL EREPORT('Emcalc_A1B',-1,cmessage)
        WRITE(6,*) '**** SRES A1B **** ',2100,' ****'
      ELSE IF (iyr1 < 198) THEN
        class(:,:) = yclass(:,:,198)
        cmessage = 'STOCHEM emissions have wrong date, using 1990'
        WRITE(6,*) cmessage
        WRITE(6,*) iyr1,iyr2,year
! DEPENDS ON: ereport
        CALL EREPORT('Emcalc_A1B',-1,cmessage)
        WRITE(6,*) '**** SRES A1B **** ',1990,' ****'
      ELSE
        class(:,:) = yclass(:,:,iyr1)+REAL(yearx-iyr1*10)*              &
     &    (yclass(:,:,iyr2)-yclass(:,:,iyr1))/10.0
        WRITE(6,*) '**** SRES A1B **** ',YEARX,' ****'
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
        IF(k /= i_o3 .AND. k /= i_hno3) emiss(k,:,:)=emiss(k,:,:)/lmolec
        IF (k == i_c5h8) THEN ! Isoprene
! DEPENDS ON: isop
          CALL ISOP(month,isopre)
          emiss(i_c5h8,:,:) = emiss(i_c5h8,:,:) * isopre(:,:)
        END IF
      END DO

      END SUBROUTINE EMCALC_A1B
#endif
