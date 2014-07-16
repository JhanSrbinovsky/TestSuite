#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!***********************************************************************
! Calculates the leaf turnover rate as a function of temperature and
! soil water availability
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!
!    Programming standard:
!
!***********************************************************************
      SUBROUTINE LEAF_LIT (LAND_PTS,VEG_PTS,VEG_INDEX,N,FSMC,TSTAR      &
     &,                    G_LEAF)

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
     & FSMC(LAND_PTS)                                                   &
                                  ! IN Soil moisture availability
!                                 !    factor.
     &,TSTAR(LAND_PTS)                                                  &
                                  ! IN Surface temperature (K).
     &,G_LEAF(LAND_PTS)                                                 &
                                  ! OUT Rate of leaf turnover
!                                 !     (/360days).
     &,FM,FT                      ! WORK Soil moisture and leaf
!                                        temperature amplifiers of
!                                        leaf turnover.

      INTEGER                                                           &
     & J,L                        ! Loop counters

#include "nstypes.h"
#include "trif.h"

!-----------------------------------------------------------------------
! Calculate the leaf turnover rate
!-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

        FT = 1.0
        FM = 1.0
        IF (TSTAR(L)  <   TLEAF_OF(N)) THEN
          FT = 1.0 + DGL_DT(N)*(TLEAF_OF(N)-TSTAR(L))
        ELSEIF (FSMC(L)  <   FSMC_OF(N)) THEN
          FM = 1.0 + DGL_DM(N)*(FSMC_OF(N)-FSMC(L))
        ENDIF

        G_LEAF(L) = G_LEAF_0(N)*FT*FM

      ENDDO

      RETURN

      END SUBROUTINE LEAF_LIT
#endif
