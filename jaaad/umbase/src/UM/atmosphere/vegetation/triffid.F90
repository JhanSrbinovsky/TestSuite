#if defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine TRIFFID ------------------------------------------------
!!!
!!!                     Top-down
!!!                     Representation of
!!!                     Interactive
!!!                     Foliage and
!!!                     Flora
!!!                     Including
!!!                     Dynamics
!!!
!!! Purpose : Simulates changes in vegetation structure, areal
!!!           coverage and the carbon contents of vegetation and soil.
!!!           can be used to advance these variables dynamically
!!!           (GAMMA=1/TIMESTEP) or to iterate towards  equilibrium
!!!           (GAMMA --> 0.0, FORW=1.0).
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4     10/97     New Deck. Peter Cox
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.
!!!                    Richard Betts
!!!  5.2   15/11/00    Re-Written for New Dynamics      M. Best
!  6.2  01/03/06  increase CS and RESP_S to dimension 4 and 5, and pass
!                 down to SOILCARB.                     C.D. Jones
!!!
!!!END ----------------------------------------------------------------
      SUBROUTINE TRIFFID (LAND_PTS,TRIF_PTS,TRIF_INDEX,FORW,GAMMA       &
     &,                   FRAC_VS,FRAC_AGRIC,G_LEAF,NPP,RESP_S          &
     &,                   RESP_W,CS,FRAC,HT,LAI,RESP_FRAC               &
     &,                   C_VEG,CV,LIT_C,LIT_C_T)



      IMPLICIT NONE

#include "nstypes.h"

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,TRIF_PTS                                                         &
                                  ! IN Number of points on which
!                                 !    TRIFFID may operate
     &,TRIF_INDEX(LAND_PTS)                                             &
                                  ! IN Indices of land points on
!                                 !    which TRIFFID may operate
     &,L,N,T                      ! WORK Loop counters

      REAL                                                              &
     & FORW                                                             &
                                  ! IN Forward timestep weighting.
     &,FRAC_VS(LAND_PTS)                                                &
                                  ! IN Total fraction of gridbox
!                                 !    covered by veg or soil.
     &,GAMMA                                                            &
                                  ! IN Inverse timestep (/360days).
     &,FRAC_AGRIC(LAND_PTS)                                             &
                                  ! IN Fraction of agriculture.
     &,G_LEAF(LAND_PTS,NPFT)                                            &
                                  ! IN Turnover rate for leaf and
!                                 !    fine root biomass (/360days).
     &,NPP(LAND_PTS,NPFT)                                               &
                                  ! INOUT Net primary productivity
!                                 !       (kg C/m2/360days).
     &,RESP_S(LAND_PTS,5)                                               &
                                 ! INOUT Soil respiration
!                                 !       (kg C/m2/360days).
     &,RESP_W(LAND_PTS,NPFT)                                            &
                                  ! INOUT Wood maintenance respiration
!                                 !       (kg C/m2/360days).
     &,CS(LAND_PTS,4)                                                   &
                                 ! INOUT Soil carbon (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                  ! INOUT Fractional cover of each
!                                 !       Functional Type.
     &,HT(LAND_PTS,NPFT)                                                &
                                  ! INOUT Vegetation height (m).
     &,LAI(LAND_PTS,NPFT)                                               &
                                  ! INOUT Leaf area index.
     &,C_VEG(LAND_PTS,NPFT)                                             &
                                  ! OUT Total carbon content of
!                                 !     the vegetation (kg C/m2).
     &,CV(LAND_PTS)                                                     &
                                  ! OUT Gridbox mean vegetation
!                                 !     carbon (kg C/m2).
     &,LIT_C(LAND_PTS,NPFT)                                             &
                                  ! OUT Carbon Litter (kg C/m2/360days).
     &,LIT_C_T(LAND_PTS)          ! OUT Gridbox mean carbon litter
!                                 !     (kg C/m2/360days).

      REAL                                                              &
     & DCVEG(LAND_PTS,NPFT)                                             &
                                  ! WORK Change in vegetation carbon
!                                 !      during the timestep
!                                 !      (kg C/m2/timestep).
     &,DFRAC(LAND_PTS,NPFT)                                             &
                                  ! WORK Change in areal fraction
!                                 !      during the timestep
!                                 !      (/timestep).
     &,FRAC_FLUX                                                        &
                                  ! WORK PFT fraction to be used
!                                 !      in the calculation of
!                                 !      the gridbox mean fluxes.
     &,RESP_FRAC(LAND_PTS)                                              &
                                  ! WORK  respired fraction of RESP_S
     &,LAI_BAL(LAND_PTS,NPFT)                                           &
                                  ! WORK Leaf area index in balanced
!                                 !      growth state.
     &,LEAF(LAND_PTS,NPFT)                                              &
                                  ! WORK Leaf biomass (kg C/m2).
     &,PC_S(LAND_PTS,NPFT)                                              &
                                  ! WORK Net carbon flux available
!                                 !      for spreading
!                                 !      (kg C/m2/yr).
     &,PHEN(LAND_PTS,NPFT)                                              &
                                  ! WORK Phenological state.
     &,ROOT(LAND_PTS,NPFT)                                              &
                                  ! WORK Root biomass (kg C/m2).
     &,WOOD(LAND_PTS,NPFT)        ! WORK Woody biomass (kg C/m2).

#include "trif.h"

!----------------------------------------------------------------------
! Loop through Functional Types
!----------------------------------------------------------------------
      DO N=1,NPFT

!----------------------------------------------------------------------
! Loop through TRIFFID points
!----------------------------------------------------------------------
        DO T=1,TRIF_PTS
          L=TRIF_INDEX(T)

!----------------------------------------------------------------------
! Diagnose the balanced-growth leaf area index and the associated leaf,
! wood, root and total vegetation carbon
!----------------------------------------------------------------------
          LAI_BAL(L,N) = (A_WS(N)*ETA_SL(N)*HT(L,N)                     &
     &              /A_WL(N))**(1.0/(B_WL(N)-1))
          LEAF(L,N) = SIGL(N)*LAI_BAL(L,N)
          ROOT(L,N) = LEAF(L,N)
          WOOD(L,N) = A_WL(N)*(LAI_BAL(L,N)**B_WL(N))
          C_VEG(L,N) = LEAF(L,N) + ROOT(L,N) + WOOD(L,N)
!----------------------------------------------------------------------
! Diagnose the phenological state
!----------------------------------------------------------------------
          PHEN(L,N) = LAI(L,N)/LAI_BAL(L,N)

        ENDDO

!----------------------------------------------------------------------
! Update vegetation carbon contents
!----------------------------------------------------------------------
! DEPENDS ON: vegcarb
        CALL VEGCARB (LAND_PTS,TRIF_PTS,TRIF_INDEX,N,FORW               &
     &,               GAMMA,G_LEAF(1,N),NPP(1,N),RESP_W(1,N)            &
     &,               LEAF(1,N),ROOT(1,N),WOOD(1,N)                     &
     &,               DCVEG(1,N),PC_S(1,N))

      ENDDO

!-----------------------------------------------------------------------
! Diagnose the new value of Canopy Height, Leaf Area Index and Total
! Vegetation Carbon
!-----------------------------------------------------------------------
      DO N=1,NPFT

        DO T=1,TRIF_PTS
          L=TRIF_INDEX(T)

          HT(L,N) = WOOD(L,N) / (A_WS(N) * ETA_SL(N))                   &
     &            * (A_WL(N)/WOOD(L,N))**(1.0/B_WL(N))
          LAI_BAL(L,N) = LEAF(L,N) / SIGL(N)
          LAI(L,N) = PHEN(L,N) * LAI_BAL(L,N)
          C_VEG(L,N) = LEAF(L,N) + ROOT(L,N) + WOOD(L,N)

        ENDDO

      ENDDO

!----------------------------------------------------------------------
! Update the areal coverage of each functional type
!----------------------------------------------------------------------
! DEPENDS ON: lotka
      CALL LOTKA (LAND_PTS,TRIF_PTS,TRIF_INDEX                          &
     &,           C_VEG,FORW,FRAC_VS,FRAC_AGRIC,GAMMA,LAI_BAL,PC_S      &
     &,           FRAC,DFRAC)

!----------------------------------------------------------------------
! Diagnose the litterfall from the carbon balance of each vegetation
! type
!----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)

        LIT_C_T(L) = 0.0

        DO N=1,NPFT
          FRAC_FLUX=FRAC(L,N)-(1.0-FORW)*DFRAC(L,N)
          LIT_C(L,N) = NPP(L,N)-GAMMA/FRAC_FLUX*(C_VEG(L,N)*FRAC(L,N)   &
     &               -(C_VEG(L,N)-DCVEG(L,N))*(FRAC(L,N)-DFRAC(L,N)))
          LIT_C_T(L) = LIT_C_T(L)+FRAC_FLUX*LIT_C(L,N)
        ENDDO
      ENDDO

!----------------------------------------------------------------------
! Call SOIL_C to update the soil carbon content
!----------------------------------------------------------------------
! DEPENDS ON: soilcarb
      CALL SOILCARB (LAND_PTS,TRIF_PTS,TRIF_INDEX                       &
     &,              LIT_C,FRAC,FRAC_AGRIC,RESP_FRAC                    &
     &,              FORW,GAMMA,LIT_C_T,RESP_S,CS)

!----------------------------------------------------------------------
! Diagnose the gridbox mean vegetation carbon
!----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)
        CV(L) = 0.0
        DO N=1,NPFT
          CV(L) = CV(L)+FRAC(L,N)*C_VEG(L,N)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE TRIFFID
#endif
