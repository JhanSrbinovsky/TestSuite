#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine CANCAP ------------------------------------------------
!!!
!!! Purpose : Calculate the heat capacity of a given PFT from its LAI
!!!
!!! version for CRAY YMP
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  5.2  15/11/00 New deck.  M. Best
!!!  5.4  11/04/02   Add can_model=4 option for canopy snow model with
!!!                  needleleaf tree tile.  R. Essery
!!!
!!!END -----------------------------------------------------------------
!**********************************************************************
      SUBROUTINE CANCAP (LAND_PTS,VEG_PTS,VEG_INDEX,CAN_MODEL,FT        &
     &,                  HT,LAI,CANHC,VFRAC)

      IMPLICIT NONE

       INTEGER                                                          &
     & LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_PTS)                                              &
                                  ! IN Index of vegetated points.
     &,CAN_MODEL                  ! IN Swith for thermal vegetation
!                                 !    canopy

      INTEGER                                                           &
     & FT                         ! IN Plant functional type.

      REAL                                                              &
     & HT(LAND_PTS)                                                     &
                                  ! IN Vegetation height (m).
     &,LAI(LAND_PTS)                                                    &
                                  ! IN Leaf area index.
     &,CANHC(LAND_PTS)                                                  &
                                  ! OUT Areal heat capacity of
!                                 !     vegetation canopy (J/K/m2).
     &,VFRAC(LAND_PTS)            ! OUT Fractional canopy coverage.

      REAL                                                              &
     & LAI_BAL(LAND_PTS)                                                &
                                  ! WORK Leaf area index in balanced
!                                 !      growth state.
     &,LEAF(LAND_PTS)                                                   &
                                  ! WORK Leaf biomass (kg C/m2).
     &,WOOD(LAND_PTS)             ! WORK Woody biomass (kg C/m2).

      INTEGER                                                           &
     & J,L                        ! WORK Loop counters.

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL                                                              &
     & HLEAF                                                            &
                                  ! Specific heat capacity of leaves
!                                 ! (J / K / kg Carbon).
     &,HWOOD                      ! Specific heat capacity of wood
!                                 ! (J / K / kg Carbon).
      PARAMETER ( HLEAF=5.7E4, HWOOD=1.1E4 )

#include "nstypes.h"
#include "pftparm.h"
#include "trif.h"

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        CANHC(L) = 0.
        VFRAC(L) = 0.
      ENDDO

      IF (CAN_MODEL  ==  2) THEN
!     Radiative canopy without heat capacity
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        CANHC(L) = 0.
        VFRAC(L) = 1. - EXP(-KEXT(FT)*LAI(L))
      ENDDO

      ELSEIF (CAN_MODEL == 3 .OR. CAN_MODEL == 4) THEN
!     Radiative canopy with heat capacity
        DO J=1,VEG_PTS
          L = VEG_INDEX(J)
          LAI_BAL(L) = ( A_WS(FT)*ETA_SL(FT)*HT(L) /                    &
     &                   A_WL(FT) )**(1.0/(B_WL(FT)-1))
          LEAF(L) = SIGL(FT)*LAI_BAL(L)
          WOOD(L) = A_WL(FT)*(LAI_BAL(L)**B_WL(FT))
          CANHC(L) = HLEAF*LEAF(L) + HWOOD*WOOD(L)
          VFRAC(L) = 1. - EXP(-KEXT(FT)*LAI(L))
        ENDDO

      ENDIF

      RETURN
      END SUBROUTINE CANCAP
#endif
