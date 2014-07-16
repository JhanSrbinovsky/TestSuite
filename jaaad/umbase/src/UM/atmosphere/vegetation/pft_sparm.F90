#if defined(A19_1A) || defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Routine to calculate the land surface parameters of a given PFT from
! its areal fraction and structural properties.
!
! Written by Peter Cox (June 1997)
!!!  5.2    15/11/00   Re-Written for New Dynamics      M. Best
!**********************************************************************
      SUBROUTINE PFT_SPARM  (LAND_PTS,N,TILE_INDEX,TILE_PTS             &
     &,                      HT,LAI,SATCON,CATCH_T,INFIL_T,Z0_T)


      IMPLICIT NONE

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                  ! IN Number of land points.
     &,N                                                                &
                                  ! IN Plant functional type.
     &,TILE_PTS                                                         &
                                  ! IN Number of land points which
!                                 !    include the surface type.
     &,TILE_INDEX(LAND_PTS)                                             &
                                  ! IN Indices of land points which
!                                 !    include the surface type.
     &,J,L                        ! WORK Loop counters.

      REAL                                                              &
     & HT(LAND_PTS)                                                     &
                                  ! IN Vegetation height (m).
     &,LAI(LAND_PTS)                                                    &
                                  ! IN Leaf area index.
     &,SATCON(LAND_PTS)                                                 &
                                  ! IN Saturated hydraulic conductivity
!                                 !    (kg/m2/s).
     &,CATCH_T(LAND_PTS)                                                &
                                  ! OUT Canopy capacity (kg/m2).
     &,INFIL_T(LAND_PTS)                                                &
                                  ! OUT Maximum surface infiltration
!                                 !     rate (kg/m2/s).
     &,Z0_T(LAND_PTS)             ! OUT Roughness length (m).

#include "nstypes.h"
#include "pftparm.h"

      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        Z0_T(L) = DZ0V_DH(N) * HT(L)
        CATCH_T(L) = CATCH0(N) + DCATCH_DLAI(N) * LAI(L)
        INFIL_T(L) = INFIL_F(N) * SATCON(L)
      ENDDO


      RETURN
      END SUBROUTINE PFT_SPARM
#endif
