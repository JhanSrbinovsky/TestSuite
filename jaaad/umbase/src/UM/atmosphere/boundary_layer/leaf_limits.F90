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
!**********************************************************************
      SUBROUTINE LEAF_LIMITS(LAND_FIELD,VEG_PTS,VEG_INDEX,FT            &
     &,                      NL,DQ,APAR,TL,CA,OA,PSTAR,FSMC             &
     &,                 CLOS_PTS,OPEN_PTS,CLOS_INDEX,OPEN_INDEX         &
     &,                 CI,RD,WCARB,WEXPT,WLITE ,FAPARV_layer,VCMAX)


      IMPLICIT NONE

#include "nstypes.h"
#include "trif.h"
#include "c_0_dg_c.h"


!  Factors in expressions for limitation of photosynthesis by transport
!  of products, for C3 and C4 respectively. Clim.Dyn. 15, 183-203 (1999)
!  Appendix A.
      real, PARAMETER      ::  fwe_c3 = 0.5
      real, PARAMETER      ::  fwe_c4 = 20000.0
      
!  q10 factor for plant respiration.
      real, PARAMETER      ::  q10_leaf = 2.0
!  (mol/sec) / (watts) conversion for PAR:
      real, PARAMETER      ::  CONPAR = 2.19e5
      
      REAL                                                              &
     & FD(NPFT)                                                         &
                              ! Dark respiration coefficient.
     & ,NEFF(NPFT)            ! Constant relating VCMAX and leaf N

      DATA FD      /   0.015,  0.015,   0.015,     0.025, 0.015 /
      DATA NEFF    /  0.8E-3 , 0.8E-3 , 0.8E-3 ,  0.4E-3, 0.8E-3  /


      INTEGER                                                           &
     & LAND_FIELD                                                       &
                                  ! IN Total number of land points.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)                                            &
                                  ! IN Index of vegetated points
                                  !    on the land grid.
     &,FT                         ! IN Plant functional type.

      INTEGER                                                           &
     & CLOS_INDEX(LAND_FIELD)                                           &
                                  ! OUT Index of land points
                                  !     with closed stomata.
     &,CLOS_PTS                                                         &
                                  ! OUT Number of land points
                                  !     with closed stomata.
     &,OPEN_INDEX(LAND_FIELD)                                           &
                                  ! OUT Index of land points
                                  !     with open stomata.
     &,OPEN_PTS                   ! OUT Number of land points
!                                 !     with open stomata.

      REAL                                                              &
     & NL(LAND_FIELD)                                                   &
                                  ! IN Leaf nitrogen
!                                 !    concentration (kg N/kg C).
     &,DQ(LAND_FIELD)                                                   &
                                  ! IN Canopy level specific humidity
!                                 !    deficit (kg H2O/kg air).
     &,APAR(LAND_FIELD)                                                 &
                                  ! IN Absorbed PAR (W/m2)
     &,TL(LAND_FIELD)                                                   &
                                  ! IN Leaf temperature (K).
     &,CA(LAND_FIELD)                                                   &
                                  ! IN Canopy CO2 pressure (Pa).
     &,OA(LAND_FIELD)                                                   &
                                  ! IN Atmospheric O2 pressure (Pa).
     &,PSTAR(LAND_FIELD)                                                &
                                  ! IN Atmospheric pressure (Pa).
     &,FSMC(LAND_FIELD)                                                 &
                                  ! IN Soil water factor.
     &,GL(LAND_FIELD)                                                   &
                                  ! OUT Leaf conductnace for H2O (m/s).
     &,AL(LAND_FIELD)                                                   &
                                  ! OUT Net Leaf photosynthesis
!                                 !     (mol CO2/m2/s).
     &,CI(LAND_FIELD)                                                   &
                                  ! OUT Internal CO2 pressure (Pa).
     &,RD(LAND_FIELD)                                                   &
                                  ! OUT Dark respiration (mol CO2/m2/s).
     &,WCARB(LAND_FIELD)                                                &
                                  ! OUT Carboxylation, ...
     &,WLITE(LAND_FIELD)                                                &
                                  !     ... Light, and ...
     &,WEXPT(LAND_FIELD)                                                &
                                  !     ... export limited gross ...
!                                 !     ... photosynthetic rates ...
!                                 !     ... (mol CO2/m2/s).
     &,ACR(LAND_FIELD)                                                  &
                                  ! WORK Absorbed PAR
!                                 !      (mol photons/m2/s).
     &,CCP(LAND_FIELD)                                                  &
                                  ! WORK Photorespiratory compensatory
!                                 !      point (mol/m3).
     &,CONV(LAND_FIELD)                                                 &
                                  ! WORK Factor for converting mol/m3
!                                 !      into Pa (J/mol).
     &,DENOM(LAND_FIELD)                                                &
                                  ! WORK Denominator in equation for VCM
     &,GLCO2(LAND_FIELD)                                                &
                                  ! WORK Leaf conductnace for CO2 (m/s).
     &,KC(LAND_FIELD)                                                   &
                                  ! WORK Michaelis constant for CO2 (Pa)
     &,KO(LAND_FIELD)                                                   &
                                  ! WORK Michaelis constant for O2 (Pa).
     &,QTENF(LAND_FIELD)                                                &
                                  ! WORK Q10 function.
     &,TAU(LAND_FIELD)                                                  &
                                  ! WORK CO2/O2 specificity ratio.
     &,TDEGC(LAND_FIELD)                                                &
                                  ! WORK Leaf temperature (deg C).
     &,VCM(LAND_FIELD)                                                  &
                                  ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco (mol CO2/m2/s).
     &,VCMAX(LAND_FIELD)                                                &
                                  ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco - without the
!                                 !      temperature factor
!                                 !      (mol CO2/m2/s).

!     LINA
     &, Vm(LAND_FIELD)                                                  &
     &, FAPARV_layer(LAND_FIELD)                                        &
                                 ! IN Profile of absorbed PAR.
     &, RES(LAND_FIELD)                                                 &
                                 ! OUT light inhibited leaf resp.(mol CO2/m2/s)
     &, Leaf_RD(LAND_FIELD)

      INTEGER                                                           &
     & J,L                        ! WORK Loop counters.


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

        IF (FSMC(L)==0.0 .OR. DQ(L) >= DQCRIT(FT)                       &
     &                     .OR. APAR(L)==0.0) THEN
          CLOS_PTS = CLOS_PTS + 1
          CLOS_INDEX(CLOS_PTS) = J
        ELSE
          OPEN_PTS = OPEN_PTS + 1
          OPEN_INDEX(OPEN_PTS) = J
        ENDIF


      ENDDO

!----------------------------------------------------------------------
! Calculate the photosynthetic parameters
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        VCMAX(L) = NEFF(FT) * NL(L)
         TDEGC(L) = TL(L) - ZeroDegC

! TAU is the Rubisco specificity for CO2 relative to O2. The numbers
! in this equation are from Cox, HCTN 24, "Description ... Vegetation
! Model", equation 53.
        TAU(L) = 2600.0 * (0.57 ** (0.1 * (TDEGC(L) - 25.0)))
        CCP(L) = 0.5 * OA(L) / TAU(L) * C3(FT)


        QTENF(L) = VCMAX(L) * (q10_leaf ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))              &
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)  ! Cox, HCTN 24, equation 49.
        !RD(L) = FD(FT) * VCM(L)
        LEAf_RD(L) = FD(FT) * VCM(L)
        !lina

      ENDDO

!DIR$ IVDEP
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        VCMAX(L) = NEFF(FT) * NL(L)
         TDEGC(L) = TL(L) - ZeroDegC

! TAU is the Rubisco specificity for CO2 relative to O2. The numbers
! in this equation are from Cox, HCTN 24, "Description ... Vegetation
! Model", equation 53.
         TAU(L) = 2600.0 * (0.57 ** (0.1 * (TDEGC(L) - 25.0)))
         CCP(L) = 0.5 * OA(L) / TAU(L) * REAL( C3(FT) )

         QTENF(L) = VCMAX(L) * (q10_leaf ** (0.1 * (TDEGC(L) - 25.0)))
         DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))             &
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)  ! Cox, HCTN 24, equation 49.

        RD(L) = FD(FT) * VCM(L)

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
        CI(L) = (CA(L) - CCP(L)) * F0(FT)                               &
     &        * (1 - DQ(L) / DQCRIT(FT)) + CCP(L)

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
        ACR(L) = APAR(L) / CONPAR

! Mercado et al, Tellus 59B, 553-565 (2007), page 560:
        IF (ACR(L)*1.E6 *FAPARV_layer(L) >  10.) then
          RD(L)=( 0.5-0.05 *log(ACR(L)*FAPARV_layer(L)*1.e6))*Leaf_RD(l)
        ELSE
          RD(L)=Leaf_RD(L)
        ENDIF

      ENDDO

!----------------------------------------------------------------------
! Calculate the gross photosynthesis for RuBP-Carboxylase, Light and
! Export limited photosynthesis (Collatz et al., 1992).
!----------------------------------------------------------------------
!DIR$ IVDEP
      IF (C3(FT)==1) THEN

        DO J=1,OPEN_PTS
          L = VEG_INDEX(OPEN_INDEX(J))
! The numbers
! in these 2 equations are from Cox, HCTN 24, "Description ... Vegetation
! Model", equations 54 and 55.
          KC(L) = 30.0 * (2.1 ** (0.1 * (TDEGC(L) - 25.0)))
          KO(L) = 30000.0 * (1.2 ** (0.1 * (TDEGC(L) - 25.0)))

          WCARB(L) = VCM(L) * (CI(L) - CCP(L))                          &
     &             / (CI(L) + KC(L) * (1. + OA(L) / KO(L)))

          WLITE(L) = ALPHA(FT) * ACR(L) * (CI(L) - CCP(L))              &
     &             / (CI(L) + 2 * CCP(L))

          WEXPT(L) = fwe_c3 * VCM(L)
        ENDDO

      ELSE

        DO J=1,OPEN_PTS

          L = VEG_INDEX(OPEN_INDEX(J))

          WCARB(L) = VCM(L)

          WLITE(L) = ALPHA(FT) * ACR(L)

          WEXPT(L) = fwe_c4 * VCM(L) * CI(L) / PSTAR(L)

        ENDDO

      ENDIF

      RETURN
      END SUBROUTINE LEAF_LIMITS
#endif
