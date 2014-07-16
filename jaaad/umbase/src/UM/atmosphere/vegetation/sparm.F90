#if defined(A19_1A) || defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Routine to calculate the gridbox mean land surface parameters from
! the areal fractions of the surface types and the structural
! properties of the plant functional types.
!
! Written by Peter Cox (June 1997)
!!!  5.2    15/11/00   Re-Written for New Dynamics      M. Best
!!!  5.4    11/04/02   Canopy snow model for needleleaf tree tile.
!!!                    R. Essery
!**********************************************************************
      SUBROUTINE SPARM (LAND_PTS,NTILES,CAN_MODEL                       &
     &,                 TILE_PTS,TILE_INDEX,FRAC,HT,LAI,SATCON          &
     &,                 CATCH_SNOW,CATCH_TILE,INFIL_TILE,Z0_TILE)

      IMPLICIT NONE

#include "nstypes.h"

      INTEGER                                                           &
     & LAND_PTS                                                         &
                             ! IN Number of land points to be processed.
     &,NTILES                                                           &
                             ! IN Number of surface tiles.
     &,CAN_MODEL                                                        &
                                    ! IN Swith for thermal vegetation
     &,TILE_PTS(NTYPE)                                                  &
                                    ! IN Number of land points which
!                                   !    include the nth surface type.
     &,TILE_INDEX(LAND_PTS,NTYPE)   ! IN Indices of land points which
!                                   !    include the nth surface type.

      REAL                                                              &
     & FRAC(LAND_PTS,NTYPE)                                             &
                                  ! IN Fractional cover of each
!                                 !    surface type.
     &,HT(LAND_PTS,NPFT)                                                &
                                  ! IN Vegetation height (m).
     &,LAI(LAND_PTS,NPFT)                                               &
                                  ! IN Leaf area index.
     &,SATCON(LAND_PTS)           ! IN Saturated hydraulic
!                                 !    conductivity (kg/m2/s).

      REAL                                                              &
     & CATCH_TILE(LAND_PTS,NTILES)                                      &
                                    ! OUT Canopy capacity for each tile
!                                   !     (kg/m2).
     &,INFIL_TILE(LAND_PTS,NTILES)                                      &
                                    ! OUT Maximum surface infiltration
!                                   !     rate for each tile (kg/m2/s).
     &,Z0_TILE(LAND_PTS,NTILES)                                         &
                                    ! OUT Roughness length for each
!                                   !     tile (m).
     &,CATCH_SNOW(LAND_PTS)         ! OUT Snow capacity for NLT tile
!                                   !     (kg/m2).
      REAL                                                              &
     & CATCH(LAND_PTS)                                                  &
                                  ! WORK GBM canopy capacity (kg/m2).
     &,CATCH_T(LAND_PTS,NTYPE)                                          &
                                  ! WORK Capacities for types.
     &,FZ0(LAND_PTS)                                                    &
                                  ! WORK Aggregation function of Z0.
     &,INFIL(LAND_PTS)                                                  &
                                  ! WORK GBM infiltration rate(kg/m2/s).
     &,INFIL_T(LAND_PTS,NTYPE)                                          &
                                  ! WORK Infiltration rates for types.
     &,Z0(LAND_PTS)                                                     &
                                  ! WORK GBM roughness length (m).
     &,Z0_T(LAND_PTS,NTYPE)       ! WORK Roughness lengths for types.

      INTEGER                                                           &
     & J,L,N                      ! WORK Loop counters

#include "pftparm.h"
#include "nvegparm.h"
#include "blend_h.h"

!----------------------------------------------------------------------
! Set parameters for vegetated surface types
!----------------------------------------------------------------------
      DO N=1,NPFT
! DEPENDS ON: pft_sparm
        CALL PFT_SPARM (LAND_PTS,N,TILE_INDEX(1,N),TILE_PTS(N)          &
     &,                 HT(1,N),LAI(1,N),SATCON                         &
     &,                 CATCH_T(1,N),INFIL_T(1,N),Z0_T(1,N))
      ENDDO

      IF (CAN_MODEL  ==  4) THEN
        N=2
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          CATCH_SNOW(L) = 4.4*LAI(L,N)
        ENDDO
      ENDIF

!----------------------------------------------------------------------
! Set parameters for non-vegetated surface types
!----------------------------------------------------------------------
      DO N=NPFT+1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          CATCH_T(L,N) = CATCH_NVG(N-NPFT)
          INFIL_T(L,N) = INFIL_NVG(N-NPFT)*SATCON(L)
          Z0_T(L,N) = Z0_NVG(N-NPFT)
        ENDDO
      ENDDO

      IF (NTILES  ==  1) THEN
!----------------------------------------------------------------------
! Form means and copy to tile arrays if required for aggregate tiles
!----------------------------------------------------------------------
        DO L=1,LAND_PTS
        CATCH(L) = 0.0
        FZ0(L) = 0.0
          INFIL(L) = 0.0
          Z0(L) = 0.0
      ENDDO

      DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          FZ0(L) = FZ0(L) + FRAC(L,N) / (LOG(LB / Z0_T(L,N)))**2
        ENDDO
      ENDDO
        DO L=1,LAND_PTS
          Z0(L) = LB * EXP(-SQRT(1. / FZ0(L)))
        ENDDO

        DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          CATCH(L) = CATCH(L) + FRAC(L,N) * CATCH_T(L,N)
            INFIL(L) = INFIL(L) + FRAC(L,N) * INFIL_T(L,N)
        ENDDO
      ENDDO

        DO L=1,LAND_PTS
!         Canopy capacity is average over non-lake surface types
          CATCH_TILE(L,1) = 0.
          IF (FRAC(L,7) <  1.)                                          &
     &      CATCH_TILE(L,1) = CATCH(L) / (1. - FRAC(L,7))
          INFIL_TILE(L,1) = INFIL(L)
          Z0_TILE(L,1) = Z0(L)
        ENDDO

      ELSE
!----------------------------------------------------------------------
! Copy surface-type arrays to tiles if separate tiles used
!----------------------------------------------------------------------
        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            CATCH_TILE(L,N) = CATCH_T(L,N)
            INFIL_TILE(L,N) = INFIL_T(L,N)
            Z0_TILE(L,N) = Z0_T(L,N)
      ENDDO
        ENDDO

      ENDIF

      RETURN
      END SUBROUTINE SPARM
#endif
