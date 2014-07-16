#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Routine to calculate the bulk stomatal resistance and the canopy
! CO2 fluxes
!
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!
!    Programming standard:
!
!**********************************************************************

!***********************************************************************
! Calculates the canopy resistance, net photosynthesis and transpiration
! by scaling-up the leaf level response using the "Big-Leaf" approach
! of Sellers et al. (1994)
!
! Written by Peter Cox (May 1995)
!***********************************************************************
      SUBROUTINE CANOPY (ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX            &
     &,                  VEG_PTS,VEG_INDEX                              &
     &,                  FT,DQC,IPAR,TSTAR,CO2C,O2,PSTAR                &
     &,                  FPAR,FSMC,LAI                                  &
     &,                  GC,ANETC,CI,RDC)

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                                  ! IN Number of points on a row
     &,ROWS                                                             &
                                  ! IN Number of rows in a theta field
     &,LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,LAND_INDEX(LAND_PTS)                                             &
                                  ! IN Index of land points on the
!                                 !    P-grid.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_PTS)        ! IN Index of vegetated points
!                                 !    on the land grid.

      INTEGER                                                           &
     & FT                         ! IN Plant functional type.

      REAL                                                              &
     & CO2C(LAND_PTS)                                                   &
                                  ! IN Canopy level CO2 concentration
!                                 !    (kg CO2/kg air).
     &,DQC(LAND_PTS)                                                    &
                                  ! IN Canopy level specific humidity
!                                 !    deficit (kg H2O/kg air).
     &,O2                                                               &
                                  ! IN Atmospheric O2 concentration
!                                 !    (kg O2/kg air).
     &,PSTAR(LAND_PTS)                                                  &
                                  ! IN Surface pressure (Pa).
     &,IPAR(ROW_LENGTH,ROWS)                                            &
                                  ! IN Incident PAR (W/m2).
     &,TSTAR(LAND_PTS)                                                  &
                                  ! IN Surface temperature (K).
     &,FPAR(LAND_PTS)                                                   &
                                  ! IN PAR absorption factor.
     &,FSMC(LAND_PTS)                                                   &
                                  ! IN Soil water factor.
     &,LAI(LAND_PTS)              ! IN Leaf area index
!                                 !    (m2 leaf/m2 ground).


      REAL                                                              &
     & ANETC(LAND_PTS)                                                  &
                                  ! OUT Net canopy photosynthesis
!                                 !     (mol CO2/m2/s).
     &,CI(LAND_PTS)                                                     &
                                  ! OUT Internal CO2 concentration
!                                 !     (mol CO2/m3).
     &,GC(LAND_PTS)                                                     &
                                  ! OUT Canopy conductance for H2O
!                                 !     (m/s).
     &,RDC(LAND_PTS)                                                    &
                                  ! OUT Canopy dark respiration
!                                 !     (mol CO2/m2/s).
     &,ANETL(LAND_PTS)                                                  &
                                  ! WORK Net leaf photosynthesis
!                                 !      (mol CO2/m2/s/LAI).
     &,APAR(LAND_PTS)                                                   &
                                  ! WORK PAR absorbed by the top leaf
!                                 !      (W/m2).
     &,CA(LAND_PTS)                                                     &
                                  ! WORK Canopy level CO2 pressure
!                                 !      (Pa).
     &,DQM(LAND_PTS)                                                    &
                                  ! WORK Canopy level humidity
!                                 !      deficit (mol H2O/m3).
     &,GL(LAND_PTS)                                                     &
                                  ! WORK Leaf conductance for H2O
!                                 !      (m/s).
     &,OA(LAND_PTS)                                                     &
                                  ! WORK Atmospheric O2 pressure
!                                 !      (Pa).
     &,RD(LAND_PTS)               ! WORK Dark respiration of top leaf
!                                 !      (mol CO2/m2/s).

      INTEGER                                                           &
     & I,J,K,L                      ! WORK Loop counters.

#include "nstypes.h"
#include "trif.h"

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL                                                              &
     & R                          ! Gas constant (J/K/mol)
      PARAMETER (R = 8.3144)
#include "ccarbon.h"

!-----------------------------------------------------------------------
! Calculate the atmospheric pressures of CO2 and O2
!-----------------------------------------------------------------------
      DO K=1,VEG_PTS
        L = VEG_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH

        CA(L) = CO2C(L) / EPCO2 * PSTAR(L)
        OA(L) = O2 / EPO2 * PSTAR(L)
        DQM(L) = DQC(L) / EPSILON * PSTAR(L) / (R * TSTAR(L))

!-----------------------------------------------------------------------
! Calculate the PAR absorbed by the top leaf
!-----------------------------------------------------------------------
        APAR(L) = (1 - OMEGA(FT)) * IPAR(I,J)

      ENDDO

!-----------------------------------------------------------------------
! Call the leaf level model for the top leaf of the C3 and C4 plants
!-----------------------------------------------------------------------

      IF ( C3(FT)  ==  1 ) THEN

! DEPENDS ON: leaf_c3
        CALL LEAF_C3 (LAND_PTS,VEG_PTS,VEG_INDEX,FT                     &
     &,               DQC,APAR,TSTAR,CA,OA,PSTAR,FSMC                   &
     &,               GL,ANETL,CI,RD)

      ELSE

! DEPENDS ON: leaf_c4
        CALL LEAF_C4 (LAND_PTS,VEG_PTS,VEG_INDEX,FT                     &
     &,               DQC,APAR,TSTAR,CA,OA,PSTAR,FSMC                   &
     &,               GL,ANETL,CI,RD)

      ENDIF

!-----------------------------------------------------------------------
! Scale-up to the canopy level
!-----------------------------------------------------------------------
      DO K=1,VEG_PTS
        L = VEG_INDEX(K)

        ANETC(L) = ANETL(L) * FPAR(L)
        GC(L) = FPAR(L) * GL(L)
        RDC(L) = RD(L) * FPAR(L)

      ENDDO

      RETURN

      END SUBROUTINE CANOPY
#endif
