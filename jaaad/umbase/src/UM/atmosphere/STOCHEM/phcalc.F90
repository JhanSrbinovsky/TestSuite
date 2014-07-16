#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE PHCALC(HP,HP1,CTHNO3,CTSO2,CTNH3,SO4,LL,M,TC,KH,KE)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : To return [H+].
!-
!-   Inputs  : HP (Trial [H+]), CTs , LL,M,TC,KHs,KEs
!-   Outputs : HP1
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.4    18/03/97  Created.  C.E. Johnson
!  5.5    19/11/03  Minor improvement to calculation of kh4_eff
!                   M.G. Sanderson
!  6.1    20/10/04  Minor tidying of code. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson

!VVV  V2.2  PHCALC 15/III/00
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!-----------------------------------------------------------------------
      REAL, DIMENSION(6), INTENT(IN):: KH,KE
      REAL, INTENT(IN) :: HP,CTHNO3,CTSO2,CTNH3,SO4,LL,M,TC
      REAL, INTENT(OUT) :: HP1
      REAL, DIMENSION(12) :: CS
      REAL :: KH4_EFF,REUS

! Effective KH from Seinfeld and Pandis p350
      KH4_EFF = KH(4) * (1.0 + (KE(2)/HP) * (1.0 + KE(3)/HP))

! 1) Obtain dissolved species concentrations in (mol/l):
!    HNO3, SO2, NH3 from equilibrium equations.


      REUS=RGC*TC*1.0E6/(NA*PSTAR) ! converts molecules/cm3 to atm.
!          HNO3(aq)+NO3
      CS(2)=CTHNO3/LL
!          CS(4) is [SO2(aq)] + [HSO3-] + [SO3--].
      CS(4)=CTSO2/(LL+PSTAR/                                            &
     &      (kh4_eff*rgc*tc*1.0e3)) ! 1.0e3 converts from mol/m3 to mol/
!          NH3+NH4 (Aq)
      CS(5)=CTNH3/LL
!          CO2   (Aq)
      CS(11)=KH(6)*XX_CO2*M*REUS !
!          NO3  - HNO3 is assumed completly dissociated.
      CS(6) = CS(2)
!          HSO3
      CS(7) = CS(4)/(1.0+HP/KE(2)+KE(3)/HP)
!          SO3
      CS(8) = CS(7)*KE(3)/HP
!          NH4
      CS(9) = CS(5)/(1+KE(6)/(HP*KE(4)))
!          SO4
      CS(10) = SO4
!          HCO3
      CS(12) = CS(11)/(1+HP/KE(5))

! 2) Calculate HP using:
!    [H+]=[NO3-]+[HSO3-]+2[SO3--]+2[SO4--]+[HCO3-]-[NH4+]

      HP1=CS(6)+CS(7)+2.0*CS(8)+2.0*CS(10)+CS(12)-CS(9)

      END SUBROUTINE PHCALC
#endif
