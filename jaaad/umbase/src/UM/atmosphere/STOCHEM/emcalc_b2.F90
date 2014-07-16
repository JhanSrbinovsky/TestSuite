#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EMCALC_B2(emiss,month,year)
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
!  6.1    06/09/04  Added emissions of toluene, o-xylene and DMS. Same
!                   for each year currently. M.G. Sanderson
!  6.2    01/03/06  Removed calculation of isoprene emission as now
!                   done interactively. Added varying emissions of
!                   aromatics.  M.G. Sanderson
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
      CHARACTER(72) :: cmessage

      class = 0.0
! Note that termites and wetland CH4 emissions are specified in
! EMUPDT and EMREAD respectively.
! Paddy and animals are assumed to be part of the anthropogenic.

!
!     -1990-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/22.8  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/393.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/282.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/1.4   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/6.0   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/3.0   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/53.6  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/6.8   ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/9.3   ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/6.2   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/7.1   ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/2.5   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,     0.0, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_toluen)= (/12.3  ,     0.0,   0.0,  0.0,  0.0 /)
      class(:,i_oxyl)=   (/9.9   ,     0.0,   0.0,  0.0,  0.0 /)
      class(:,i_so2)=    (/65.1  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/70.9  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,199)=class(:,:)*1e9

!
!     -2000-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/25.4  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/536.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/307.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/1.5   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/6.8   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/3.4   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/61.1  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/7.7   ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/10.6  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/7.1   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/8.0   ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/2.8   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/69.0  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/12.4  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/10.1  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,200)=class(:,:)*1e9

!
!     -2010-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/32.5  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/635.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/354.0 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/1.7   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/7.6   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/3.8   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/68.5  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/8.7   ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/11.8  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/7.9   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/9.0   ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/3.1   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/65.9  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/14.6  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/11.4  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,201)=class(:,:)*1e9

!
!     -2020-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/36.3  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/711.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/392.5 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.5   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/10.1   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/6.0   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/75.7  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/12.7  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/16.6  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/9.1   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/13.3  ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/3.7   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/61.3  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/15.9  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/12.9  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,202)=class(:,:)*1e9

!
!     -2030-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/41.3  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/675.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/445.9 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.7   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/10.8  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/6.6   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/79.7  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/13.8  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/18.0  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/9.6   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/14.4  ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/4.0   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/60.3  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/17.6  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/14.2  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,203)=class(:,:)*1e9

!
!     -2040-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/45.7  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/768.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/473.3 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.0   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/11.8  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/7.4   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/85.2  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/15.3  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/19.8  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/10.3  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/16.0  ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/4.3   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/59.0  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/19.0  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/15.3  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,204)=class(:,:)*1e9

!
!     -2050-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/46.6  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/851.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/482.9 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.0   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/12.0  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/7.5   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/86.0  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/15.5  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/20.1  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/10.4  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/16.2  ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/4.3   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/55.7  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/19.2  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/15.5  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,205)=class(:,:)*1e9

!
!     -2060-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/48.3  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/966.0 ,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/489.4 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/3.0   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/11.8  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/7.3   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/84.9  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/15.2  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/19.7  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/10.2  ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/15.9  ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/4.3   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/53.8  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/19.0  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/15.3  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,206)=class(:,:)*1e9

!
!     -2070-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/48.5  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/1125.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/486.5 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.7   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/10.9  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/6.7   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/80.1  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/13.9  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/18.1  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/9.6   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/14.5  ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/4.0   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/50.9  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/17.9  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/14.4  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,207)=class(:,:)*1e9

!
!     -2080-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/51.4  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/1303.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/474.4 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.5   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/10.1  ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/6.0   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/75.7  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/12.7  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/16.6  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/9.0   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/13.3  ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/3.7   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/50.0  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/17.0  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/13.7  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,208)=class(:,:)*1e9

!
!     -2090-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/53.0  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/1448.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/453.5 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.2   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/9.0   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/5.3   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/70.2  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/11.2  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/14.8  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/8.3   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/11.7  ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/3.4   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/49.0  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/15.8  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/12.7  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,209)=class(:,:)*1e9

!
!     -2100-B2-----------EMISSION CLASSES-(Tg/yr)----------------------
!
!                         Anth   ,Biomass ,  Veg , Soil,Oceans
      class(:,i_no)=     (/53.3  ,     7.1,   0.0,  5.6,  2.6 /) ! NO
      class(:,i_co)=     (/1567.0,   500.0, 150.0,  0.0, 50.0 /) ! CO
      class(:,i_ch4)=    (/453.4 ,    55.0,   0.0,  0.0, 13.0 /) ! CH4
      class(:,i_hcho)=   (/2.1   ,     1.5,   0.0,  0.0,  0.0 /) ! HCHO
      class(:,i_h2)=     (/20.0  ,   20.00,   0.0,  5.0,  5.0 /) ! H2
      class(:,i_c2h6)=   (/8.5   ,     4.4,   3.5,  0.0,  0.0 /) ! C2H6
      class(:,i_ch3cho)= (/4.8   ,     4.8,   0.0,  0.0,  0.0 /) ! CH3CH
      class(:,i_nc4h10)= (/67.0  ,     2.6,   8.0,  0.0,  0.0 /) ! C4H10
      class(:,i_c2h4)=   (/10.4  ,     8.4,  20.0,  0.0,  0.0 /) ! C2H4
      class(:,i_c3h6)=   (/13.7  ,     9.5,  20.0,  0.0,  0.0 /) ! C3H6
      class(:,i_c3h8)=   (/7.9   ,     1.1,   3.5,  0.0,  0.5 /) ! C3H8
      class(:,i_ch3oh)=  (/10.8  ,     8.8,   0.0,  0.0,  0.0 /) ! METHA
      class(:,i_acetone)=(/3.2   ,     0.7,  20.0,  0.0,  0.0 /) ! ACETO
      class(:,i_c5h8)=   (/0.0   ,    0.00, 450.0,  0.0,  0.0 /) ! C5H8
      class(:,i_nh3)=    (/39.4  ,     3.5,   0.0,  2.5,  8.2 /) ! NH3
      class(:,i_so2)=    (/47.9  ,     2.2,   0.0,  0.0,  0.0 /) ! SO2
      class(:,i_toluen)= (/15.1  ,     0.0,   0.0,  0.0,  0.0 /) ! TOLUE
      class(:,i_oxyl)=   (/12.1  ,     0.0,   0.0,  0.0,  0.0 /) ! o-XYL
      class(:,i_dms)=    (/0.0   ,     0.0,   0.0,  1.0, 15.0 /) ! DMS
      class(:,i_rn222)=  (/0.0   ,     0.0,   0.0,15.E-9, 0.0 /) ! RN-22
      yclass(:,:,210)=class(:,:)*1e9

      yclass(:,:,211)=yclass(:,:,210)
      yclass(:,:,198)=yclass(:,:,199)

! Mod to add 425 Tg to anthropogenic emissions for specified year,month
!      IF(YEAR == 1990 .AND. MONTH == 1) THEN
!        yclass(1,I_CO,199) = yclass(1,I_CO,199) + 425.0*1.0E12
!        WRITE(6,*) 'YEAR/MONTH ',year,'/',month,' 500 Tg ',
!     &    'added to anthropogenic EXTRA'
!      ELSE   ! to put back to original
!        yclass(1,I_CO,199) = 393.0*1.0E12
!      ENDIF

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
        class(:,:)=yclass(:,:,210)
        cmessage = 'STOCHEM emissions have wrong date, using 2100'
        WRITE(6,*) cmessage
        WRITE(6,*) iyr1,iyr2,yearx
! DEPENDS ON: ereport
        CALL EREPORT('Emcalc_B2',-1,cmessage)
        WRITE(6,*) '**** SRES B2 **** ',2100,' ****'
      ELSE IF (iyr1 < 198) THEN
        class(:,:)=yclass(:,:,198)
        cmessage = 'STOCHEM emissions have wrong date, using 1990'
        WRITE(6,*) cmessage
        WRITE(6,*) iyr1,iyr2,yearx
! DEPENDS ON: ereport
        CALL EREPORT('Emcalc_B2',-1,cmessage)
        WRITE(6,*) '**** SRES B2 **** ',1990,' ****'
      ELSE
        class(:,:)=yclass(:,:,iyr1)+REAL(yearx-iyr1*10)*                &
     &           (yclass(:,:,iyr2)-yclass(:,:,iyr1))/10.0
        WRITE(6,*) '**** SRES B2 **** ',yearx,' ****'
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

      DO k = 1, nc
! Don't scale O3 and HNO3 (will do in STRATCALC)
        IF (k/=i_o3.AND.k/=i_hno3) emiss(k,:,:)=emiss(k,:,:)/lmolec
      END DO

      END SUBROUTINE EMCALC_B2
#endif
