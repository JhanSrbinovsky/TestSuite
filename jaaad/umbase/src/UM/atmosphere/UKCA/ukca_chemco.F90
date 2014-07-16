#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  To calculate the reaction rate co-efficients for use in the
!  Backward Euler solver
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method: Calculates rate coefficients
!                 Units : cm, molecule, s
!
!   Inputs  : tc,m,h2o,o2
!   Outputs : rc
!
! Current Code Owner:       Colin Johnson/Fiona O'Connor
!
!  Code Description:
!   Language:  FORTRAN 90 (formatted)
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      SUBROUTINE UKCA_CHEMCO(nr,n_pnts,tc,m,h2o,o2,                     &
                             clw,fcloud,fdiss,k_dms,rc)

      USE UKCA_CONSTANTS,   ONLY: avc
      USE ASAD_MOD,         ONLY: peps
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"

      INTEGER, INTENT(IN) :: nr                ! No of reactions
      INTEGER, INTENT(IN) :: n_pnts            ! No of points

      REAL, INTENT(IN)    :: tc(n_pnts)        ! Temperature
      REAL, INTENT(IN)    :: m(n_pnts)         ! Air density
      REAL, INTENT(IN)    :: h2o(n_pnts)       ! Water vapour
      REAL, INTENT(IN)    :: o2(n_pnts)        ! Oxygen
      REAL, INTENT(IN)    :: clw(n_pnts)       ! Cloud Liquid water
      REAL, INTENT(IN)    :: fcloud(n_pnts)    ! Cloud fraction

!     Dissolved fraction - allows for 2 ions:
      REAL, INTENT(IN)    :: fdiss(n_pnts,jpspec,2)

!     Rate coeffs to allow parameterisation of DMS products:
      REAL, INTENT(OUT)   :: k_dms(n_pnts,5)

!     Thermal rate coefficients:
      REAL, INTENT(OUT)   :: rc(n_pnts,nr)

!     Local variables

      INTEGER :: i, j

      REAL :: z1(n_pnts)
      REAL :: z2(n_pnts) 
      REAL :: z3(n_pnts) 
      REAL :: z4(n_pnts) 
      REAL :: ratiob2total(n_pnts) 
      REAL :: ratioa2b(n_pnts) 
      REAL :: at(7)
      REAL :: t300(n_pnts)
      REAL :: zo(n_pnts), zi(n_pnts), arg2(n_pnts)
      REAL :: zfc(n_pnts), zr(n_pnts)
      REAL :: Hplus(n_pnts)                ! Hydrogen ion concentration
      REAL :: faq(n_pnts,jpspec)           ! total dissolved fraction


      IF (L_ukca_aerchem) THEN
        faq(:,:) = fdiss(:,:,1) + fdiss(:,:,2)
      ENDIF

!     Initialise reaction rate coefficients 

      DO j = 1, nr
        DO i = 1, n_pnts 
          rc(i,j) = 0.0
        END DO
      END DO      

!     Calculate bimolecular reaction rate coefficients 
         
!     R1: HO2 + NO = OH + NO2              3.60E-12  0.00   -270.0   
!     IUPAC no change 
         
      RC(:,1) = 3.60E-12*exp(270.0/TC(:)) 
         
!     R2: HO2 + NO3 = OH + NO2             4.00E-12  0.00      0.0   
!     IUPAC no change 
         
      RC(:,2) = 4.00E-12 
         
!     R3: HO2 + O3 = OH + O2               2.03E-16  4.57   -693.0   
!     IUPAC no change 
         
      RC(:,3) = 2.03E-16*(tc(:)/300.0)**4.57*exp(693.0/tc(:)) 
         
!     R4: HO2 + HO2 = H2O2                 2.20E-13  0.00   -600.0    
!     IUPAC no change - Note 4  
          
      RC(:,4) = 2.20E-13*exp(600.0/tc(:)) 
      RC(:,4) = RC(:,4)*(1.0 +                                       & 
     &          1.4E-21*h2o(:)*exp(2200.0/tc(:)) ) 
         
!     R5: HO2 + MeOO = MeOOH               3.80E-13  0.00   -780.0 
!     IUPAC 2005 update - Note 5 - Brnch A 
         
      ratiob2total(:) = 1.0 / (1.0 + 498.0*EXP(-1160.0/tc(:)))   
      RC(:,5)         = 3.80E-13*exp(780.0/tc(:))      
      RC(:,5)         = RC(:,5)*(1.0-ratiob2total(:))   
         
!     R6: HO2 + MeOO = HCHO                3.80E-13  0.00   -780.0   
!     IUPAC 2005 update - Note 6 - Brnch B 
         
      RC(:,6) = 3.80E-13*exp(780.0/tc(:)) 
      RC(:,6) = RC(:,6)*ratiob2total(:) 
         
!     R7: HO2 + EtOO = EtOOH               3.80E-13  0.00   -900.0   
!     IUPAC no change 
         
      RC(:,7) = 3.80E-13*exp(900.0/tc(:)) 
         
!     R8: HO2 + MeCO3 = MeCO3H             2.08E-13  0.00   -980.0    
!     IUPAC 2005 update - Note 8 
         
      RC(:,8) = 2.08E-13*exp(980.0/tc(:)) 
         
!     R9: HO2 + MeCO3 = MeCO2H + O3        1.04E-13  0.00   -980.0   
!     IUPAC 2005 update - note 9  
         
      RC(:,9) = 1.04E-13*exp(980.0/tc(:)) 
         
!     R10 HO2 + MeCO3 = OH + MeOO          2.08E-13  0.00   -980.0    
!     IUPAC 2005 update - Note 10 
         
      RC(:,10) = 2.08E-13*exp(980.0/tc(:)) 
         
!     R11 HO2 + n-PrOO = n-PrOOH           1.51E-13  0.00  -1300.0   
!     MCM no change 
         
      RC(:,11) = 1.51E-13*exp(1300.0/tc(:)) 
        
!     R12 HO2 + i-PrOO = i-PrOOH           1.51E-13  0.00  -1300.0          
!     MCM no change 
         
      RC(:,12) = 1.51E-13*exp(1300.0/tc(:)) 
         
!     R13 HO2 + EtCO3 = O2 + EtCO3H        3.05E-13  0.00  -1040.0   
!     MCM no change 
               
      RC(:,13) = 3.05E-13*exp(1040.0/tc(:)) 
         
!     R14 HO2 + EtCO3 = O3 + EtCO2H        1.25E-13  0.00  -1040.0   
!     MCM no change 
         
      RC(:,14) = 1.25E-13*exp(1040.0/tc(:)) 
         
!     R15 HO2 + MeCOCH2OO = MeCOCH2OOH     1.36E-13  0.00  -1250.0    
!     MCM no change - Note 15 
         
      RC(:,15) = 1.36E-13*exp(1250.0/tc(:)) 
         
!     R16 MeOO + NO = HO2 + HCHO + NO2     2.95E-12  0.00   -285.0   
!     IUPAC no change 

      RC(:,16) = 2.95E-12*exp(285.0/tc(:)) 
         
!     R17 MeOO + NO = MeONO2               2.95E-15  0.00   -285.0   
!     IUPAC no change - note 17 
         
      RC(:,17) = 2.95E-15*exp(285.0/tc(:)) 
         
!     R18 MeOO + NO3 = HO2 + HCHO + NO2    1.30E-12  0.00      0.0   
!     IUPAC no change 
         
      RC(:,18) = 1.30E-12 
         
!     R19 MeOO + MeOO = MeOH + HCHO        1.03E-13  0.00   -365.0   
!     IUPAC no change - note 19 - Brnch A 
         
      ratiob2total(:) = 1.0/(1.0+(EXP(1300.0/tc(:)))/33.0) 
      RC(:,19)        = 1.03E-13*EXP(365.0/tc(:)) 
      RC(:,19)        = RC(:,19)*(1.0-ratiob2total(:)) 
               
!     R20 MeOO + MeOO = 2HO2 + 2HCHO       1.03E-13  0.00   -365.0   
!     IUPAC no change - note 20 - Brnch B 
         
      RC(:,20) = 1.03E-13*EXP(365.0/tc(:)) 
      RC(:,20) = RC(:,20)*ratiob2total(:) 
         
!     R21 MeOO + MeCO3 = HO2 + HCHO + MeOO 1.80E-12  0.00   -500.0   
!     IUPAC no change - note 21 
         
      RC(:,21) = 1.80E-12*EXP(500.0/tc(:)) 
         
!     R22 MeOO + MeCO3 = MeCO2H + HCHO     2.00E-13  0.00   -500.0   
!     IUPAC no change - note 22 
         
      RC(:,22) = 2.00E-13*EXP(500.0/tc(:)) 
         
!     R23 EtOO + NO = MeCHO + HO2 + NO2    2.60E-12  0.00   -380.0   
!     IUPAC no change - note 23 
         
      RC(:,23) = 2.60E-12*EXP(380.0/tc(:)) 
         
!     R24 EtOO + NO3 = MeCHO + HO2 + NO2   2.30E-12  0.00      0.0   
!     IUPAC no change 
         
      RC(:,24) = 2.30E-12 

!     R25 EtOO + MeCO3 = MeCHO + HO2 + MeOO 4.40E-13  0.00  -1070.0   
!     IUPAC no change - note 25 
         
      RC(:,25) = 4.40E-13*EXP(1070.0/tc(:)) 
         
!     R26 MeCO3 + NO = MeOO + CO2 + NO2    7.50E-12  0.00   -290.0   
!     IUPAC no change 
         
      RC(:,26) = 7.50E-12*EXP(290.0/tc(:)) 
         
!     R27 MeCO3 + NO3 = MeOO + CO2 + NO2   4.00E-12  0.00      0.0   
!     MCM no change 
         
      RC(:,27) = 4.00E-12 
         
!     R28 n-PrOO + NO = EtCHO + HO2 + NO2  2.90E-12  0.00   -350.0   
!     IUPAC no change - note 28 
         
      RC(:,28) = 2.90E-12*EXP(350.0/tc(:)) 
         
!     R29 n-PrOO + NO3 = EtCHO + HO2 + NO2 2.50E-12  0.00      0.0   
!     MCM no change 
         
      RC(:,29) = 2.50E-12 
         
!     R30 i-PrOO + NO = Me2CO + HO2 + NO2  2.70E-12  0.00   -360.0   
!     IUPAC no change - note 30 
         
      RC(:,30) = 2.70E-12*EXP(360.0/tc(:)) 
         
!     R31 i-PrOO + NO3 = Me2CO + HO2 + NO2 2.50E-12  0.00      0.0   
!     MCM no change 
         
      RC(:,31) = 2.50E-12 
         
!     R32 EtCO3 + NO = EtOO + CO2 + NO2    6.70E-12  0.00   -340.0   
!     IUPAC no change 
         
      RC(:,32) = 6.70E-12*EXP(340.0/tc(:)) 
 
!     R33 EtCO3 + NO3 = EtOO + CO2 + NO2   4.00E-12  0.00      0.0   
!     MCM no change 
         
      RC(:,33) = 4.00E-12 

!     R34 MeCOCH2OO + NO = MeCO3+HCHO+NO2  2.80E-12  0.00   -300.0   
!     Tyndall et al. - note 34 
         
      RC(:,34) = 2.80E-12*EXP(300.0/tc(:)) 
         
!     R35 MeCOCH2OO + NO3 = MeCO3+HCHO+NO2 2.50E-12  0.00      0.0   
!     MCM no change 
         
      RC(:,35) = 2.50E-12 
         
!     R36 NO + NO3 = NO2 + NO2             1.80E-11  0.00   -110.0   
!     IUPAC no change 
         
      RC(:,36) = 1.80E-11*EXP(110.0/tc(:)) 
         
!     R37 NO + O3 = NO2                    1.40E-12  0.00   1310.0   
!     IUPAC no change 
         
      RC(:,37) = 1.40E-12*EXP(-1310.0/tc(:)) 
         
!     R38 NO2 + O3 = NO3                   1.40E-13  0.00   2470.0   
!     IUPAC no change 
         
      RC(:,38) = 1.40E-13*EXP(-2470.0/tc(:)) 
         
!     R39 NO3 + HCHO = HONO2 + HO2 + CO    2.00E-12  0.00   2440.0   
!     IUPAC no change - note 39 
         
      RC(:,39) = 2.00E-12*EXP(-2440.0/tc(:)) 
         
!     R40 NO3 + MeCHO = HONO2 + MeCO3      1.40E-12  0.00   1860.0   
!     IUPAC no change 
         
      RC(:,40) = 1.40E-12*EXP(-1860.0/tc(:)) 
         
!     R41 NO3 + EtCHO = HONO2 + EtCO3      3.46E-12  0.00   1862.0   
!     MCM no change - note 41 
         
      RC(:,41) = 3.46E-12*EXP(-1862.0/tc(:))     

!     R42 NO3 + Me2CO = HONO2 + MeCOCH2OO  3.00E-17  0.00      0.0   
!     IUPAC no change - note 42 
         
      RC(:,42) = 3.00E-17 
        
!     R43 N2O5 + H2O = HONO2 + HONO2       2.50E-22  0.00      0.0   
!     IUPAC no change - note 43 
         
      RC(:,43) = 2.50E-22 
         
!     R44 O(3P) + O3 = O2 + O2             8.00E-12  0.00   2060.0   
!     IUPAC no change 
         
      RC(:,44) = 8.00E-12*EXP(-2060.0/tc(:)) 
         
!     R45 O(1D) + CH4 = OH + MeOO          1.05E-10  0.00      0.0   
!     IUPAC no change - note 45 

      RC(:,45) = 1.05E-10 
         
!     R46 O(1D) + CH4 = HCHO + H2          7.50E-12  0.00      0.0   
!     IUPAC no change - note 46 
         
      RC(:,46) = 7.50E-12 
         
!     R47 O(1D) + CH4 = HCHO + HO2 + HO2   3.45E-11  0.00      0.0   
!     IUPAC no change - note 47 
         
      RC(:,47) = 3.45E-11 
         
!     R48 O(1D) + H2O = OH + OH            2.20E-10  0.00      0.0   
!     IUPAC no change 
         
      RC(:,48) = 2.20E-10 
         
!     R49 O(1D) + N2 = O(3P) + N2          2.10E-11  0.00   -115.0   
!     Ravishankara et al. - note 49 
         
      RC(:,49) = 2.10E-11*EXP(115.0/tc(:)) 
         
!     R50 O(1D) + O2 = O(3P) + O2          3.20E-11  0.00    -67.0   
!     IUPAC no change 
         
      RC(:,50) = 3.20E-11*EXP(67.0/tc(:)) 
         
!     R51 OH + CH4 = H2O + MeOO            1.85E-12  0.00   1690.0   
!     IUPAC no change 
         
      RC(:,51) = 1.85E-12*EXP(-1690.0/tc(:)) 
        
!     R52 OH + C2H6 = H2O + EtOO           6.90E-12  0.00   1000.0   
!     IUPAC no change 
         
      RC(:,52) = 6.90E-12*EXP(-1000.0/tc(:)) 
         
!     R53 OH + C3H8 = n-PrOO + H2O         7.60E-12  0.00    585.0   
!     IUPAC no change - note 53 - Brnch A 
         
      ratioa2b(:) = 226.0 * tc(:)**(-0.64)*EXP(-816.0/tc(:)) 
      RC(:,53) = 7.60E-12*EXP(-585.0/tc(:)) 
      RC(:,53) = RC(:,53)*(ratioa2b(:)/(ratioa2b(:)+1.0)) 
               
!     R54 OH + C3H8 = i-PrOO + H2O         7.60E-12  0.00    585.0   
!     IUPAC no change - note 54 - Brnch B 
         
      RC(:,54) = 7.60E-12*EXP(-585.0/tc(:)) 
      RC(:,54) = RC(:,54)/(ratioa2b(:)+1.0) 
               
!     R55 OH + CO = HO2                    1.44E-13  0.00      0.0   
!     IUPAC 2005 update - note 55 
         
      RC(:,55) = 1.44E-13*(1.0+m(:)/4.2E19) 
         
!     R56 OH + EtCHO = H2O + EtCO3         5.10E-12  0.00   -405.0   
!     IUPAC no change  
         
      RC(:,56) = 5.10E-12*EXP(405.0/tc(:)) 
         
!     R57 OH + EtOOH = H2O + MeCHO + OH    8.01E-12  0.00      0.0   
!     MCM no change 
         
      RC(:,57) = 8.01E-12 
         
!     R58 OH + EtOOH = H2O + EtOO          1.90E-12  0.00   -190.0   
!     MCM no change 
         
      RC(:,58) = 1.90E-12*EXP(190.0/tc(:)) 
         
!     R59 OH + H2 = H2O + HO2              7.70E-12  0.00   2100.0   
!     IUPAC no change 
         
      RC(:,59) = 7.70E-12*EXP(-2100.0/tc(:))    
      
!     R60 OH + H2O2 = H2O + HO2            2.90E-12  0.00    160.0  
!     IUPAC no change

      RC(:,60) = 2.90E-12*EXP(-160.0/tc(:)) 

      IF (L_ukca_aerchem) THEN

!       Reduce rate to account for dissolved H2O2
        RC(:,60) = RC(:,60)*(1.0-(faq(:,13)*fcloud(:)))

      ENDIF

!     R61 OH + HCHO = H2O + HO2 + CO       5.40E-12  0.00   -135.0  
!     IUPAC 2004 update - note 61

      RC(:,61) = 5.40E-12*EXP(135.0/tc(:))

!     R62 OH + HO2 = H2O                   4.80E-11  0.00   -250.0  
!     IUPAC no change

      RC(:,62) = 4.80E-11*EXP(250.0/tc(:))

!     R63 OH + HO2NO2 = H2O + NO2          1.90E-12  0.00   -270.0  
!     IUPAC no change

      RC(:,63) = 1.90E-12*EXP(270.0/tc(:))

!     R64 OH + HONO2 = H2O + NO3           1.50E-13  0.00      0.0  
!     IUPAC no change - note 64

      z1(:)    = 2.4e-14 * exp( 460.0/tc(:)) 
      z3(:)    = 6.5e-34 * exp(1335.0/tc(:)) 
      z4(:)    = 2.7e-17 * exp(2199.0/tc(:)) 
      z2(:)    = z3(:)*m(:) / (1.0+z3(:)*m(:)/z4(:)) 
      RC(:,64) = z1(:) + z2(:) 

!     R65 OH + HONO = H2O + NO2            2.50E-12  0.00   -260.0  
!     IUPAC no change

      RC(:,65) = 2.50E-12*exp(260.0/tc(:)) 

!     R66 OH + MeOOH = H2O + HCHO + OH     1.02E-12  0.00   -190.0  
!     IUPAC no change - note 66

      RC(:,66) = 1.02E-12*EXP(190.0/tc(:)) 

!     R67 OH + MeOOH = H2O + MeOO          1.89E-12  0.00   -190.0  
!     IUPAC no change - note 67

      RC(:,67) = 1.89E-12*EXP(190.0/tc(:)) 

!     R68 OH + MeONO2 = HCHO + NO2 + H2O   4.00E-13  0.00    845.0  
!     IUPAC no change

      RC(:,68) = 4.00E-13*EXP(-845.0/tc(:)) 

!     R69 OH + Me2CO = H2O + MeCOCH2OO     8.80E-12  0.00   1320.0  
!     IUPAC no change - note 69

      RC(:,69) = 8.80E-12*EXP(-1320.0/tc(:)) 
      
!     R70 OH + Me2CO = H2O + MeCOCH2OO     1.70E-14  0.00   -420.0  
!     IUPAC no change - note 70

      RC(:,70) = 1.70E-14*EXP(420.0/tc(:))     

!     R71 OH + MeCOCH2OOH = H2O+MeCOCH2OO  1.90E-12  0.00   -190.0  
!     MCM no change - note 71

      RC(:,71) = 1.90E-12*EXP(190.0/tc(:)) 

!     R72 OH + MeCOCH2OOH = OH + MGLY      8.39E-12  0.00      0.0  
!     MCM no change - note 72

      RC(:,72) = 8.39E-12 
      
!     R73 OH + MeCHO = H2O + MeCO3         4.40E-12  0.00   -365.0  
!     IUPAC no change

      RC(:,73) = 4.40E-12*EXP(365.0/tc(:)) 

!     R74 OH + NO3 = HO2 + NO2             2.00E-11  0.00      0.0  
!     IUPAC no change

      RC(:,74) = 2.00E-11

!     R75 OH + O3 = HO2 + O2               1.70E-12  0.00    940.0  
!     IUPAC no change

      RC(:,75) = 1.70E-12*EXP(-940.0/tc(:))

!     R76 OH + OH = H2O + O(3P)            6.31E-14  2.60   -945.0  
!     IUPAC no change - note 76

      RC(:,76) = 6.31E-14*((tc(:)/300.0)**2.6)*EXP(945.0/tc(:))

!     R77 OH + PAN = HCHO + NO2 + H2O      3.00E-14  0.00      0.0  
!     IUPAC no change - note 77

      RC(:,77) = 3.00E-14

!     R78 OH + PPAN = MeCHO + NO2 + H2O    1.27E-12  0.00      0.0  
!     MCM no change

      RC(:,78) = 1.27E-12

!     R79 OH + n-PrOOH = n-PrOO + H2O      1.90E-12  0.00   -190.0  
!     MCM no change

      RC(:,79) = 1.90E-12*EXP(190.0/tc(:))

!     R80 OH + n-PrOOH = EtCHO + H2O + OH  1.10E-11  0.00      0.0  
!     MCM no change

      RC(:,80) = 1.10E-11

!     R81 OH + i-PrOOH = i-PrOO + H2O      1.90E-12  0.00   -190.0  
!     MCM no change

      RC(:,81) = 1.90E-12*EXP(190.0/tc(:))

!     R82 OH + i-PrOOH = Me2CO + OH        1.66E-11  0.00      0.0  
!     MCM no change

      RC(:,82) = 1.66E-11

!     R83 O(3P) + NO2 = NO + O2            5.50E-12  0.00   -188.0  
!     IUPAC no change     

      RC(:,83) = 5.50E-12*EXP(188.0/tc(:))

!     R84 HO2S + O3S = HO2S + O2           2.03E-16  4.57   -693.0  
!     IUPAC no change - note 84

      RC(:,84) = 2.03E-16*((tc(:)/300.0)**4.57)*EXP(693.0/tc(:))

!     R85 OHS + O3S = OHS + O2             1.70E-12  0.00    940.0  
!     IUPAC no change

      RC(:,85) = 1.70E-12*EXP(-940.0/tc(:))

!     R86 O(1D)S + H2O = H2O               2.20E-10  0.00      0.0  
!     IUPAC no change

      RC(:,86) = 2.20E-10

!     R87 O(1D)S + N2 = O(3P)S + N2        2.10E-11  0.00   -115.0  
!     Ravishankara et al. - note 49

      RC(:,87) = 2.10E-11*EXP(115.0/tc(:))

!     R88 O(1D)S + O2 = O(3P)S + O2        3.20E-11  0.00    -67.0  
!     IUPAC no change

      RC(:,88) = 3.20E-11*EXP(67.0/tc(:))
      
!     R89 HO2 + HO2 = H2O2 + O2     0.00 1.90E-33  0.00   -980.0 0.00E+00  0.00      0.0  
!     IUPAC no change - note 1

      at(1) = 0.00
      at(2) = 1.90E-33
      at(3) = 0.00
      at(4) = -980.0
      at(5) = 0.00E+00
      at(6) = 0.00
      at(7) = 0.0

      t300(:) = tc(:)/300.0
      zo(:)   = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
      zi(:)   = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))
              
      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,89) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,89) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((log10(zr(i)))**2))
          RC(i,89) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO      
      RC(:,89) = RC(:,89)*(1.0 + 1.4E-21*h2o(:)*exp(2200.0/tc(:)) )

      
!     R90 HO2 + NO2 = HO2NO2 + m    0.60 1.80E-31 -3.20  0.0 4.70E-12  0.00      0.0  
!     IUPAC no change       

      at(1) = 0.60
      at(2) = 1.80E-31
      at(3) = -3.20
      at(4) = 0.0
      at(5) = 4.70E-12
      at(6) = 0.00
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))
      
      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,90) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,90) = zo(i)
        ELSE
          IF( at(1) <= 1.0 )THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = EXP( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,90) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R91 HO2NO2 + m = HO2 + NO2  0.60 4.10E-05  0.00  10650.0 4.80E+15  0.00  11170.0  
!     IUPAC no change      

      at(1) = 0.60
      at(2) = 4.10E-5
      at(3) = 0.0
      at(4) = 10650.0
      at(5) = 4.80E+15
      at(6) = 0.00
      at(7) = 11170.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))
      
      DO i=1,n_pnts
        IF( zo(i) < peps )THEN
          RC(i,91) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,91) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,91) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R92 MeCO3 + NO2 = PAN + m  0.30 2.70E-28 -7.10   0.0 1.20E-11 -0.90      0.0  
!     IUPAC no change      

      at(1) = 0.30
      at(2) = 2.70E-28
      at(3) = -7.10
      at(4) = 0.0
      at(5) = 1.21E-11
      at(6) = -0.90
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))
      
      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,92) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,92) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,92) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R93 PAN + m = MeCO3 + NO2    0.30 4.90E-03  0.00  12100.0 5.40E+16  0.00  13830.0  
!     IUPAC no change

      at(1) = 0.30
      at(2) = 4.90E-03
      at(3) = 0.00
      at(4) = 12100.0
      at(5) = 5.40E+16
      at(6) = 0.00
      at(7) = 13830.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))
      
      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,93) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,93) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,93) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R94 N2O5 + m = NO2 + NO3  0.35 1.30E-03 -3.50  11000.0 9.70E+14  0.10  11080.0  
!     IUPAC no change - note 6

      at(1) = 0.35
      at(2) = 1.30E-03
      at(3) = -3.50
      at(4) = 11000.0
      at(5) = 9.70E+14
      at(6) = 0.10
      at(7) = 11080.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,94) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,94) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,94) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R95 NO2 + NO3 = N2O5 + m   0.35 3.60E-30 -4.10   0.0 1.90E-12  0.20      0.0  
!     IUPAC no change - note 7    

      at(1) = 0.35
      at(2) = 3.60E-30
      at(3) = -4.10
      at(4) = 0.0
      at(5) = 1.90E-12
      at(6) = 0.20
      at(7) = 0.0
      
      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,95) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,95) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,95) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R96 O(3P) + O2 = O3 + m      0.00 5.70E-34 -2.60  0.0 0.00E+00  0.00      0.0  
!     IUPAC no change - note 8       

      at(1) = 0.00
      at(2) = 5.70E-34
      at(3) = -2.60
      at(4) = 0.0
      at(5) = 0.00E+00
      at(6) = 0.00
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF ( zo(i) < peps ) THEN
          RC(i,96) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,96) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,96) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

! R97 OH + NO = HONO + m       1420.00 7.40E-31 -2.40      0.0 3.30E-11 -0.30      0.0  
! IUPAC no change - note 9   

      at(1) = 1420.0
      at(2) = 7.40E-31
      at(3) = -2.40
      at(4) = 0.0
      at(5) = 3.30E-11
      at(6) = -0.30
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,97) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,97) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,97) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R98 OH + NO2 = HONO2 + m   0.40 3.30E-30 -3.00   0.0 4.10E-11  0.00      0.0  
!     IUPAC no change - note 10 

      at(1) = 0.40
      at(2) = 3.30E-30
      at(3) = -3.00
      at(4) = 0.0
      at(5) = 4.10E-11
      at(6) = 0.00
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,98) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,98) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,98) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R99 OH + OH = H2O2 + m      0.50 6.90E-31 -0.80   0.0 2.60E-11  0.00      0.0  
!     IUPAC no change

      at(1) = 0.50
      at(2) = 6.90E-31
      at(3) = -0.80
      at(4) = 0.0
      at(5) = 2.60E-11
      at(6) = 0.00
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))
      
      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,99) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,99) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,99) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R100 EtCO3 + NO2 = PPAN + m    0.30 2.70E-28 -7.10   0.0 1.20E-11 -0.90      0.0  
!     MCM v3.1 no change - note 12   

      at(1) = 0.30
      at(2) = 2.70E-28
      at(3) = -7.10
      at(4) = 0.0
      at(5) = 1.20E-11
      at(6) = -0.90
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,100) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,100) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)     = zo(i) / zi(i)
          arg2(i)   = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,100) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R101 PPAN + m = EtCO3 + NO2  0.36 1.70E-03  0.00  11280.0 8.30E+16  0.00  13940.0  
!     IUPAC no change  

      at(1) = 0.36
      at(2) = 1.70E-03
      at(3) = 0.00
      at(4) = 11280.0
      at(5) = 8.30E+16
      at(6) = 0.00
      at(7) = 13940.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,101) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,101) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)     = zo(i) / zi(i)
          arg2(i)   = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,101) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R102 O(3P)S + O2 = O3S + m   0.00 5.70E-34 -2.60   0.0 0.00E+00  0.00      0.0  
!     IUPAC no change - note 8 

      at(1) = 0.00
      at(2) = 5.70E-34
      at(3) = -2.60
      at(4) = 0.0
      at(5) = 0.00E+00
      at(6) = 0.00
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF ( zo(i) < peps ) THEN
          RC(i,102) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,102) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)     = zo(i) / zi(i)
          arg2(i)   = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,102) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

      IF(L_ukca_aerchem) THEN

!       R103  SO2      OH     m    H2SO4   Rate data from NASA/JPL.

        at(1) = 0.60
        at(2) = 3.00E-31
        at(3) = -3.30
        at(4) = 0.0
        at(5) = 1.50E-12
        at(6) = 0.00
        at(7) = 0.0

        zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
        zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

        DO i=1,n_pnts
          IF(zo(i) < peps ) THEN
             RC(i,103) = zi(i)
          ELSE IF( zi(i) < peps ) THEN
            RC(i,103) = zo(i)
          ELSE
            IF( at(1).le.1.0 ) THEN
               zfc(i) = at(1)
            ELSE
              zfc(i) = EXP( -tc(i)/at(1) )
            END IF
            zr(i) = zo(i) / zi(i)
            arg2(i) = 1.0/(1.0+((ALOG10(zr(i)))**2))
            rc(i,103) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*ALOG(zfc(i)))
            rc(i,103) = RC(i,103) * (1.0 - (faq(i,48) * fcloud(i)))
          END IF
        ENDDO

!       R104:  HSO3(-)(aq) + H2O2(aq) =>  H2SO4
!       Rate data from Eqn. 7 of Bower et al (1991, Atm.Env. 25A, 2401-2418)
!       Adjustments to allow for use of aqueous chemistry in model with
!       gas-phase reactions from Bergel et al (2004, JGR 109).

        Hplus(:) = 1.0E-05
        DO i=1, n_pnts
          IF(clw(i) < 1E-10) THEN
            rc(i,104) = 0.0
          ELSE
            rc(i,104) = 1.10E+12 * (EXP(-3655.3 / tc(i)))               &
                               * Hplus(i) / (0.1 + Hplus(i))
            rc(i,104) = rc(i,104) * fcloud(i) * fdiss(i,48,2)           &
                     * fdiss(i,13,1) * 1000.0 / (avc * clw(i))
          ENDIF
        ENDDO

!       DMS Oxidation from IUPAC March 2007.
!       R105 : DMS + OH = DMSO2

        rc(:,105) = 1.12E-11*EXP(-250.0/TC(:))

!       R106 : DMS + OH = DMSO2

        rc(:,106) = 9.3E-39*o2(:)*EXP(5270.0/tc(:))/                  &
                  (1.0 + 7.4E-29*o2(:)*EXP(5610.0/tc(:)))

!       R107 : DMS + NO3 = DMSO2 + HNO3

        rc(:,107) = 1.9E-13*EXP(520.0/tc(:))

!       These are not used directly, but to calculate product ratios
!       CH3SO2 => CH3 + SO2

        k_dms(:,1) = 2.6E11*EXP(-9056.0/tc(:))     ! Ayers et al 1996

!       CH3SO2 + O3 => CH3SO3

        k_dms(:,2) = 5.0E-15                       ! Yin et al 1990

!       CH3SO2 + NO2 => CH3SO3

        k_dms(:,3) = 2.2E-12                       ! Ray et al 1996

!       CH3SO3 + HO2 => MSA

        k_dms(:,4) = 4.0E-11                       ! Koga & Tanaka 1999

!       CH3SO3 => CH3 + SO3

        k_dms(:,5) = 1.1e17*EXP(-12057.0/tc(:))    ! Ayers et al 1996

      END IF    ! L_ukca_aerchem


      RETURN
      END SUBROUTINE UKCA_CHEMCO
#endif
