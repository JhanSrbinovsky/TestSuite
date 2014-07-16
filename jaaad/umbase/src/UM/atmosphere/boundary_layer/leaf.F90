#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Subroutine to calculate gridbox mean values of surface conductance
! and carbon fluxes. Also returns net primary productivity, leaf
! turnover and wood respiration of each plant functional type for
! driving TRIFFID.
!
!
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1992) C3 photosynthesis model, and the
!      Collatz et al. (1991) C4 photosynthesis model.
! (ii) Jacobs (1994) CI/CA closure.
! Written by Peter Cox (February 1996)
! Adapted for MOSES II tile model by Richard Essery (July 1997)
! Coded into UMVN 6.2  by Pete Falloon (July 2006)
!**********************************************************************
      SUBROUTINE LEAF (LAND_FIELD,VEG_PTS,VEG_INDEX,FT                  &
     &,                CLOS_PTS,OPEN_PTS,CLOS_INDEX,OPEN_INDEX          &
     &,                FSMC,TL,CA,CI,RD,WCARB,WEXPT,WLITE               &
     &,                GL,AL)
     
      IMPLICIT NONE

#include "nstypes.h"
#include "trif.h"

      INTEGER                                                           &
     & LAND_FIELD                                                       &
                                  ! IN Total number of land points.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)                                            &
                                  ! IN Index of vegetated points
                                  !    on the land grid.
     &,FT                                                               &
                                  ! IN Plant functional type.
     &,CLOS_INDEX(LAND_FIELD)                                           &
                                  ! IN Index of land points
                                  !    with closed stomata.
     &,CLOS_PTS                                                         &
                                  ! IN Number of land points
                                  !    with closed stomata.
     &,OPEN_INDEX(LAND_FIELD)                                           &
                                  ! IN Index of land points
                                  !    with open stomata.
     &,OPEN_PTS                   ! IN Number of land points
                                  !    with open stomata.

      REAL                                                              &
     & FSMC(LAND_FIELD)                                                 &
                                  ! IN Soil water factor.
     &,TL(LAND_FIELD)                                                   &
                                  ! IN Leaf temperature (K).
     &,RD(LAND_FIELD)                                                   &
                                  ! IN Dark respiration (mol CO2/m2/s).
     &,CA(LAND_FIELD)                                                   &
                                  ! IN Canopy CO2 pressure (Pa).
     &,CI(LAND_FIELD)                                                   &
                                  ! IN Internal CO2 pressure (Pa).
     &,WCARB(LAND_FIELD)                                                &      
     &,WLITE(LAND_FIELD)                                                &
     &,WEXPT(LAND_FIELD)                                                &
                                    ! IN Carboxylation,
                                    !    Light, and
                                    !    export limited gross
                                    !    photosynthetic rates
                                    !    (mol CO2/m2/s).

     &,GL(LAND_FIELD)                                                   &
                                  ! OUT Leaf conductance for H2O (m/s).
     &,AL(LAND_FIELD)                                                   &
                                  ! OUT Net Leaf photosynthesis
!                                 !     (mol CO2/m2/s).
     &,B1(LAND_FIELD)                                                   &
     &,B2(LAND_FIELD)                                                   &
     &,B3(LAND_FIELD)                                                   &
                                  ! WORK Coefficients of the quadratic.
     &,CONV(LAND_FIELD)                                                 &
                                  ! WORK Factor for converting mol/m3
!                                 !      into Pa (J/mol).
     &,GLCO2(LAND_FIELD)                                                &
                                  ! WORK Leaf conductance for CO2 (m/s).
     &,WL(LAND_FIELD)                                                   &
                                  ! WORK Gross leaf phtosynthesis
!                                 !      (mol CO2/m2/s).
     &,WP(LAND_FIELD)             ! WORK Smoothed minimum of
!                                 !      Carboxylation and Light
!                                 !      limited gross photosynthesis
!                                 !      (mol CO2/m2/s).


      INTEGER                                                           &
     & J,L                        ! WORK Loop counters.

!-------------------------------------------------------------------
! Parameters
!-------------------------------------------------------------------
      REAL                                                              &
     & BETA1, BETA2                                                     &
                      ! Coupling coefficients for co-limitation.
     &,R                                                                &
                      ! Gas constant (J/K/mol).
     &,RATIO          ! Ratio of leaf resistance for CO2 to leaf
                     ! resistance for H2O.
      PARAMETER (BETA1 = 0.83, BETA2 = 0.93                             &
     &,          R = 8.3144, RATIO = 1.6)

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
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(L)

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

        IF (GL(L).LE.GLMIN(FT) .OR. AL(L).LE.0.0) THEN
          GL(L) = GLMIN(FT)
          GLCO2(L) = GL(L) / RATIO
          AL(L) = -RD(L)*FSMC(L)
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
        AL(L) = -RD(L) *FSMC(L)
      ENDDO

      RETURN
      END SUBROUTINE LEAF

#endif
