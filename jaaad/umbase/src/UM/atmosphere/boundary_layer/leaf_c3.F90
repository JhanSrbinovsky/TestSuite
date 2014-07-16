#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1992) C3 photosynthesis model
! (ii) Jacobs (1994) CI/CA closure.
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!
!    Programming standard:
!
!**********************************************************************
      SUBROUTINE LEAF_C3 (LAND_PTS,VEG_PTS,VEG_INDEX,FT                 &
     &,                   DQ,APAR,TL,CA,OA,PSTAR,FSMC                   &
     &,                   GL,AL,CI,RD)

      IMPLICIT NONE

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_PTS)                                              &
                                  ! IN Index of vegetated points
!                                 !    on the land grid.
     &,FT                         ! IN Plant functional type.

      REAL                                                              &
     & DQ(LAND_PTS)                                                     &
                                  ! IN Canopy level specific humidity
!                                 !    deficit (kg H2O/kg air).
     &,APAR(LAND_PTS)                                                   &
                                  ! IN Absorbed PAR (W/m2)
     &,TL(LAND_PTS)                                                     &
                                  ! IN Leaf temperature (K).
     &,CA(LAND_PTS)                                                     &
                                  ! IN Canopy CO2 pressure (Pa).
     &,OA(LAND_PTS)                                                     &
                                  ! IN Atmospheric O2 pressure (Pa).
     &,PSTAR(LAND_PTS)                                                  &
                                  ! IN Atmospheric pressure (Pa).
     &,FSMC(LAND_PTS)                                                   &
                                  ! IN Soil water factor.
     &,GL(LAND_PTS)                                                     &
                                  ! OUT Leaf conductnace for H2O (m/s).
     &,AL(LAND_PTS)                                                     &
                                  ! OUT Net Leaf photosynthesis
!                                 !     (mol CO2/m2/s).
     &,RD(LAND_PTS)                                                     &
                                  ! OUT Dark respiration (mol CO2/m2/s).
     &,CI(LAND_PTS)                                                     &
                                  ! OUT Internal CO2 pressure (Pa).
     &,ACR(LAND_PTS)                                                    &
                                  ! WORK Absorbed PAR
!                                 !      (mol photons/m2/s).
     &,B1(LAND_PTS)                                                     &
                                  !
     &,B2(LAND_PTS)                                                     &
                                  !
     &,B3(LAND_PTS)                                                     &
                                  ! WORK Coefficients of the quadratic.
     &,CCP(LAND_PTS)                                                    &
                                  ! WORK Photorespiratory compensatory
!                                 !      point (mol/m3).
     &,CONV(LAND_PTS)                                                   &
                                  ! WORK Factor for converting mol/m3
!                                 !      into Pa (J/mol).
     &,DENOM(LAND_PTS)                                                  &
                                  ! WORK Denominator in equation for VCM
     &,GLCO2(LAND_PTS)                                                  &
                                  ! WORK Leaf conductnace for CO2 (m/s).
     &,KC(LAND_PTS)                                                     &
                                  ! WORK Michaelis constant for CO2 (Pa)
     &,KO(LAND_PTS)                                                     &
                                  ! WORK Michaelis constant for O2 (Pa).
     &,QTENF(LAND_PTS)                                                  &
                                  ! WORK Q10 function.
     &,TAU(LAND_PTS)                                                    &
                                  ! WORK CO2/O2 specificity ratio.
     &,TDEGC(LAND_PTS)                                                  &
                                  ! WORK Leaf temperature (deg C).
     &,VCM(LAND_PTS)                                                    &
                                  ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco (mol CO2/m2/s).
     &,VCMAX(LAND_PTS)                                                  &
                                  ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco - without the
!                                 !      temperature factor
!                                 !      (mol CO2/m2/s).
     &,WL(LAND_PTS)                                                     &
                                  ! WORK Gross leaf phtosynthesis
!                                 !      (mol CO2/m2/s).
     &,WCARB(LAND_PTS)                                                  &
                                  ! WORK Carboxylation,
     &,WLITE(LAND_PTS)                                                  &
                                  !      Light, and
     &,WEXPT(LAND_PTS)                                                  &
                                  !      export limited gross
!                                 !      photosynthetic rates
!                                 !      (mol CO2/m2/s).
     &,WP(LAND_PTS)               ! WORK Smoothed minimum of
!                                 !      Carboxylation and Light
!                                 !      limited gross photosynthesis
!                                 !      (mol CO2/m2/s).

      INTEGER                                                           &
     & CLOS_INDEX(LAND_PTS)                                             &
                                  ! WORK Index of land points
!                                 !      with closed stomata.
     &,CLOS_PTS                                                         &
                                  ! WORK Number of land points
!                                 !      with closed stomata.
     &,OPEN_INDEX(LAND_PTS)                                             &
                                  ! WORK Index of land points
!                                 !      with open stomata.
     &,OPEN_PTS                   ! WORK Number of land points
!                                 !      with open stomata.

      INTEGER                                                           &
     & J,L                        ! WORK Loop counters.

#include "nstypes.h"
#include "trif.h"

!-------------------------------------------------------------------
! Parameters
!-------------------------------------------------------------------
      REAL                                                              &
     & BETA1, BETA2                                                     &
                      ! Coupling coefficients for co-limitation.
     &,FDC3, FDC4                                                       &
                      ! Dark respiration coefficients for C3, C4
     &,NEFFC3, NEFFC4                                                   &
                      ! Constant relating VCMAX and leaf N
!                     ! from Schulze et al. 1994 (AMAX = 0.4E-3 * NL
!                     ! - assuming dry matter is 40% carbon by mass)
!                     ! and Jacobs 1994:
!                     ! C3 : VCMAX = 2 * AMAX ; C4 : VCMAX = AMAX
!                     ! (mol/m2/s)
     &,R                                                                &
                      ! Gas constant (J/K/mol).
     &,RATIO                                                            &
                      ! Ratio of leaf resistance for CO2 to leaf
!                     ! resistance for H2O.
     &,ZERODEGC       ! Zero Celsius (K).

      PARAMETER (BETA1 = 0.83, BETA2 = 0.93                             &
     &,          FDC3 = 0.015,  FDC4 = 0.025                            &
     &,          NEFFC3 = 0.8E-3, NEFFC4 = 0.4E-3                       &
     &,          R = 8.3144  , RATIO = 1.6                              &
     &,          ZERODEGC = 273.15)

!----------------------------------------------------------------------
! Initialise counters
!----------------------------------------------------------------------
      CLOS_PTS = 0
      OPEN_PTS = 0

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

!----------------------------------------------------------------------
! Calculate the points with closed stomata
!----------------------------------------------------------------------
        IF (FSMC(L) == 0.0 .OR. DQ(L) >= DQCRIT(FT)                     &
     &                     .OR. APAR(L) == 0.0) THEN
          CLOS_PTS = CLOS_PTS + 1
          CLOS_INDEX(CLOS_PTS) = J
        ELSE
          OPEN_PTS = OPEN_PTS + 1
          OPEN_INDEX(OPEN_PTS) = J
        ENDIF

!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(L)

      ENDDO

!----------------------------------------------------------------------
! Calculate the photosynthetic parameters
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        VCMAX(L) = NEFFC3 * NL0(FT)
        TDEGC(L) = TL(L) - ZERODEGC

        TAU(L) = 2600.0 * (0.57 ** (0.1 * (TDEGC(L) - 25.0)))
        CCP(L) = 0.5 * OA(L) / TAU(L)


        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))              &
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)
        RD(L) = FDC3 * VCM(L)

      ENDDO

!DIR$ IVDEP
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        VCMAX(L) = NEFFC3 * NL0(FT)
        TDEGC(L) = TL(L) - ZERODEGC

        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))              &
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)
        RD(L) = FDC3 * VCM(L)

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
        CI(L) = (CA(L) - CCP(L)) * F0(FT)                               &
     &        * (1 - DQ(L) / DQCRIT(FT)) + CCP(L)

!----------------------------------------------------------------------
! Convert absorbed PAR into mol PAR photons/m2/s
!----------------------------------------------------------------------
        ACR(L) = APAR(L) / 2.19E5

      ENDDO

!----------------------------------------------------------------------
! Calculate the gross photosynthesis for RuBP-Carboxylase, Light and
! Export limited photosynthesis (Collatz et al., 1992).
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        KC(L) = 30.0 * (2.1 ** (0.1 * (TDEGC(L) - 25.0)))
        KO(L) = 30000.0 * (1.2 ** (0.1 * (TDEGC(L) - 25.0)))

        WCARB(L) = VCM(L) * (CI(L) - CCP(L))                            &
     &           / (CI(L) + KC(L) * (1. + OA(L) / KO(L)))

        WLITE(L) = ALPHA(FT) * ACR(L) * (CI(L) - CCP(L))                &
     &           / (CI(L) + 2 * CCP(L))

        WEXPT(L) = 0.5 * VCM(L)

      ENDDO

!----------------------------------------------------------------------
! Calculate the co-limited rate of gross photosynthesis
!----------------------------------------------------------------------

!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA1
        B2(L) = - (WCARB(L) + WLITE(L))
        B3(L) = WCARB(L) * WLITE(L)

        WP(L) = -B2(L)/(2*B1(L))                                        &
     &         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA2
        B2(L) = - (WP(L) + WEXPT(L))
        B3(L) = WP(L) * WEXPT(L)

        WL(L) = -B2(L)/(2*B1(L))                                        &
     &         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the net rate of photosynthesis
!----------------------------------------------------------------------
        AL(L) = (WL(L) - RD(L)) * FSMC(L)

!----------------------------------------------------------------------
! Diagnose the leaf conductance
!----------------------------------------------------------------------
        GLCO2(L) = (AL(L) * CONV(L)) / (CA(L) - CI(L))
        GL(L) = RATIO * GLCO2(L)

      ENDDO

!----------------------------------------------------------------------
! Close stomata at points with negative or zero net photosynthesis
! or where the leaf resistance exceeds its maximum value.
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        IF (GL(L) <= GLMIN(FT) .OR. AL(L) <= 0.0) THEN
          GL(L) = GLMIN(FT)
          GLCO2(L) = GL(L) / RATIO
          AL(L) = -RD(L) * FSMC(L)
          CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)
        ENDIF

      ENDDO

!----------------------------------------------------------------------
! Define fluxes and conductances for points with closed stomata
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        GL(L) = GLMIN(FT)
        GLCO2(L) = GL(L) / RATIO
        AL(L) = -RD(L) * FSMC(L)
        CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)

      ENDDO

      RETURN
      END SUBROUTINE LEAF_C3

!**********************************************************************
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1991) C4 photosynthesis model
! (ii) Jacobs (1994) CI/CA closure.
!
! Written by Peter Cox (February 1996)
!**********************************************************************
#endif
