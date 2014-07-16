#if defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine GROWTH -------------------------------------------------
!!!
!!! Purpose : Increments leaf, root and wood carbon.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4     10/97     New Deck. Peter Cox
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.
!!!                    Richard Betts
!!!  5.2    15/11/00   Re-Written for New Dynamics      M. Best
!!!
!!!END ----------------------------------------------------------------
      SUBROUTINE GROWTH (LAND_PTS,TRIF_PTS,TRIF_INDEX                   &
     &,                  N,DPCG_DLAI,FORW,GAMMA,PC_G                    &

     &,                  LEAF,ROOT,WOOD)

      IMPLICIT NONE

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,TRIF_PTS                                                         &
                                  ! IN Number of points on which
!                                 !    TRIFFID may operate
     &,TRIF_INDEX(LAND_PTS)                                             &
                                  ! IN Indices of land points on
!                                 !    which TRIFFID may operate
     &,N                                                                &
                                  ! IN Vegetation type.
     &,L,T                        ! WORK Loop counters

      REAL                                                              &
     & DPCG_DLAI(LAND_PTS)                                              &
                                  ! IN Rate of change of PC_G with
!                                 !    leaf area index
!                                 !    (kg C/m2/360days/LAI).
     &,FORW                                                             &
                                  ! IN Forward timestep weighting.
     &,GAMMA                                                            &
                                  ! IN Inverse timestep (/360days).
     &,PC_G(LAND_PTS)                                                   &
                                  ! IN Net carbon flux available
!                                 !    for growth (kg C/m2/360days).
     &,LEAF(LAND_PTS)                                                   &
                                  ! INOUT Leaf biomass (kg C/m2).
     &,ROOT(LAND_PTS)                                                   &
                                  ! INOUT Root biomass (kg C/m2).
     &,WOOD(LAND_PTS)             ! INOUT Woody biomass (kg C/m2).

      REAL                                                              &
     & DENOM                                                            &
                                  ! WORK Denominator of update
!                                 !      equation.
     &,DLEAF,DROOT,DWOOD                                                &
                                  ! WORK Increments to leaf, root
!                                 !      and woody biomass (kg C/m2).
     &,DL_DW                                                            &
                                  ! WORK Rate of change of leaf
!                                 !      carbon with wood carbon.
     &,DLAI_DW                                                          &
                                  ! WORK Rate of change of leaf area
!                                 !      index with wood carbon
!                                 !      (LAI m2/kg C).
     &,DR_DW                                                            &
                                  ! WORK Rate of change of root
!                                 !      carbon with wood carbon.
     &,NUMER                                                            &
                                  ! WORK Numerator of the update
!                                 !      equation.
     &,WOOD_MAX                                                         &
                                  ! WORK Maximum wood carbon (kg C/m2).
     &,WOOD_MIN                   ! WORK Minimum wood carbon (kg C/m2).

#include "nstypes.h"
#include "trif.h"
#include "descent.h"

      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)

!----------------------------------------------------------------------
! Calculate the increment to the wood carbon
!----------------------------------------------------------------------
        DL_DW = LEAF(L)/(B_WL(N)*WOOD(L))
        DR_DW = DL_DW
        DLAI_DW = DL_DW/SIGL(N)

        NUMER = PC_G(L)
        DENOM = (1+DL_DW+DR_DW)*GAMMA-FORW*DLAI_DW*DPCG_DLAI(L)
        DENOM = MAX(DENOM,DENOM_MIN)

        DWOOD = NUMER/DENOM

!----------------------------------------------------------------------
! Ensure that the local leaf area index does not drop below its
! minimum value or exceed its maximum value.
!----------------------------------------------------------------------
        WOOD_MIN = A_WL(N)*LAI_MIN(N)**B_WL(N)
        WOOD_MAX = A_WL(N)*LAI_MAX(N)**B_WL(N)
        DWOOD = MAX((WOOD_MIN-WOOD(L)),DWOOD)
        DWOOD = MIN((WOOD_MAX-WOOD(L)),DWOOD)

!----------------------------------------------------------------------
! Diagnose the increments to leaf and root carbon
!----------------------------------------------------------------------
        DLEAF = SIGL(N)*((WOOD(L)+DWOOD)/A_WL(N))**(1.0/B_WL(N))        &
     &         -LEAF(L)
        DROOT = DLEAF

!----------------------------------------------------------------------
! Update carbon contents
!----------------------------------------------------------------------
        LEAF(L) = LEAF(L)+DLEAF
        ROOT(L) = ROOT(L)+DROOT
        WOOD(L) = WOOD(L)+DWOOD

      ENDDO

      RETURN
      END SUBROUTINE GROWTH
#endif
