#if defined(A19_1A) || defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine PHENOL -------------------------------------------------
!!!
!!!
!!! Purpose :  Parametrizes leaf phenological changes and updates the
!!!            leaf area index and the leaf turnover rate.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4     10/97     New Deck. Peter Cox
!!!  5.2    15/11/00   Re-Written for New Dynamics      M. Best
!!!
!!!END -----------------------------------------------------------------
      SUBROUTINE PHENOL (LAND_PTS,VEG_PTS,VEG_INDEX,N,G_LEAF,HT         &
     &,                  DTIME_PHEN,G_LEAF_PHEN,LAI)

      IMPLICIT NONE

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_PTS)                                              &
                                  ! IN Index of vegetated points
!                                 !    on the land grid.
     &,N                          ! IN Plant functional type.

      REAL                                                              &
     & G_LEAF(LAND_PTS)                                                 &
                                  ! IN Rate of leaf turnover (/360days).
     &,HT(LAND_PTS)                                                     &
                                  ! IN Canopy height (m).
     &,DTIME_PHEN                                                       &
                                  ! IN Timestep (years).
     &,G_LEAF_PHEN(LAND_PTS)                                            &
                                  ! OUT Rate of leaf turnover
!                                 !     including leaf phenology
!                                 !     (/360days).
     &,LAI(LAND_PTS)                                                    &
                                  ! INOUT Leaf area index.
     &,DPHEN                                                            &
                                  ! WORK Increment to phenological
!                                 !      state.
     &,LAI_BAL(LAND_PTS)                                                &
                                  ! WORK Balanced growth LAI.
     &,PHEN(LAND_PTS)             ! WORK Phenological state.

      INTEGER                                                           &
     & J,L                        ! Loop counters

#include "nstypes.h"
#include "trif.h"

!-----------------------------------------------------------------------
! Diagnose the phenological state
!-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        LAI_BAL(L) = (A_WS(N)*ETA_SL(N)*HT(L)                           &
     &               /A_WL(N))**(1.0/(B_WL(N)-1))
        PHEN(L) = LAI(L)/LAI_BAL(L)
      ENDDO

!-----------------------------------------------------------------------
! Update the phenological state and output the leaf turnover rate in
! terms of the balanced growth LAI
!-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

        IF (G_LEAF(L) >  2*G_LEAF_0(N)) THEN
          DPHEN = -DTIME_PHEN*G_GROW(N)
          DPHEN = MAX(DPHEN,(0.01-PHEN(L)))
          G_LEAF_PHEN(L) = -DPHEN/DTIME_PHEN
        ELSE
          DPHEN = DTIME_PHEN*G_GROW(N)*(1.0-PHEN(L))
          DPHEN = MIN(DPHEN,(1.0-PHEN(L)))
          G_LEAF_PHEN(L) = PHEN(L)*G_LEAF(L)
        ENDIF

!-----------------------------------------------------------------------
! Update the leaf area index
!-----------------------------------------------------------------------
        PHEN(L) = PHEN(L) + DPHEN
        LAI(L) = PHEN(L)*LAI_BAL(L)

      ENDDO

      RETURN

      END SUBROUTINE PHENOL
#endif
