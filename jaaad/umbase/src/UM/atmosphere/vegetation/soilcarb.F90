#if defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine SOILCARB -----------------------------------------------
!!!
!!! Purpose : Updates carbon contents of the soil.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4     10/97     New Deck. Peter Cox
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.
!!!                    Richard Betts
!!!  5.2   15/11/00    Re-Written for New Dynamics      M. Best
!!!  6.2   01/03/06    Upgrade to include RothC.        Chris Jones.
!!!
!!!END ----------------------------------------------------------------
      SUBROUTINE SOILCARB (LAND_PTS,TRIF_PTS,TRIF_INDEX                 &
     &,                    LIT_C,FRAC,FRAC_AGRIC,RESP_FRAC              &
     &,                    FORW,GAMMA,LIT_C_T,RESP_S,CS)

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
     &,L,T                        ! WORK Loop counters

      REAL                                                              &
     & FORW                                                             &
                                  ! IN Forward timestep weighting.
     &,GAMMA                                                            &
                                  ! IN Inverse timestep (/360days).
     &,LIT_C_T(LAND_PTS)                                                &
                                  ! IN Total carbon litter
!                                 !    (kg C/m2/360days).
     &,LIT_C(LAND_PTS,NPFT)                                             &
                                  ! IN Carbon Litter (kg C/m2/360days).
     &,RESP_S(LAND_PTS,5)                                               &
                                  ! INOUT Soil respiration
!                                 !    (kg C/m2/360days).
     &,CS(LAND_PTS,4)                                                   &
                                  ! INOUT Soil carbon (kg C/m2).
     &,DCS(LAND_PTS,4)                                                  &
                                  ! WORK Increment to the soil carbon
!                                 !      (kg C/m2).
     &,DPC_DCS(LAND_PTS,4)                                              &
                                  ! WORK Rate of change of PC with
!                                 !      soil carbon (/360days).
     &,PC(LAND_PTS,4)                                                   &
                                  ! WORK Net carbon accumulation in
!                                 !      the soil (kg C/m2/360days).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                  ! INOUT Fractional cover of each
     &,FRAC_AGRIC(LAND_PTS)                                             &
                                  ! INOUT Fractional cover of each
     &,DPM_RATIO(LAND_PTS)                                              &
                                  ! WORK ratio of litter input into DPM
     &,RESP_FRAC(LAND_PTS)        ! respired fraction of RESP_S


! sum total respiration
      DO L=1,LAND_PTS
        RESP_S(L,5) = RESP_S(L,1) + RESP_S(L,2) +                       &
     &                RESP_S(L,3) + RESP_S(L,4)
      ENDDO

! calculate DPM:RPM ratio of input litter Carbon
! DEPENDS ON: dpm_rpm
      CALL DPM_RPM(LAND_PTS,TRIF_PTS,TRIF_INDEX,FRAC,                   &
     &             FRAC_AGRIC,LIT_C,DPM_RATIO)

      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)

!----------------------------------------------------------------------
! Diagnose the net local carbon flux into the soil
!----------------------------------------------------------------------
        PC(L,1) = DPM_ratio(L)*LIT_C_T(L) - RESP_S(L,1)
        PC(L,2) = (1-DPM_ratio(L))*LIT_C_T(L) - RESP_S(L,2)
        PC(L,3) = 0.46*RESP_FRAC(L)*RESP_S(L,5) - RESP_S(L,3)
        PC(L,4) = 0.54*RESP_FRAC(L)*RESP_S(L,5) - RESP_S(L,4)

!----------------------------------------------------------------------
! Variables required for the implicit and equilibrium calculations
!----------------------------------------------------------------------
        DPC_DCS(L,1) = RESP_S(L,1)/CS(L,1)
        DPC_DCS(L,2) = RESP_S(L,2)/CS(L,2)
        DPC_DCS(L,3) = RESP_S(L,3)/CS(L,3)
        DPC_DCS(L,4) = RESP_S(L,4)/CS(L,4)

!----------------------------------------------------------------------
! Save current value of soil carbon
!----------------------------------------------------------------------
        DCS(L,1) = CS(L,1)
        DCS(L,2) = CS(L,2)
        DCS(L,3) = CS(L,3)
        DCS(L,4) = CS(L,4)

      ENDDO


!----------------------------------------------------------------------
! Update soil carbon
!----------------------------------------------------------------------
! DEPENDS ON: decay
      CALL DECAY (LAND_PTS,TRIF_PTS,TRIF_INDEX                          &
     &,           DPC_DCS,FORW,GAMMA,PC,CS)

!----------------------------------------------------------------------
! Apply implicit correction to the soil respiration rate.
!----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)

        DCS(L,1) = CS(L,1) - DCS(L,1)
        DCS(L,2) = CS(L,2) - DCS(L,2)
        DCS(L,3) = CS(L,3) - DCS(L,3)
        DCS(L,4) = CS(L,4) - DCS(L,4)

        RESP_S(L,1) = RESP_S(L,1) + FORW*DPC_DCS(L,1)*DCS(L,1)
        RESP_S(L,2) = RESP_S(L,2) + FORW*DPC_DCS(L,2)*DCS(L,2)
        RESP_S(L,3) = RESP_S(L,3) + FORW*DPC_DCS(L,3)*DCS(L,3)
        RESP_S(L,4) = RESP_S(L,4) + FORW*DPC_DCS(L,4)*DCS(L,4)

      ENDDO

! sum total respiration
      DO L=1,LAND_PTS
        RESP_S(L,5) = RESP_S(L,1) + RESP_S(L,2) +                       &
     &                RESP_S(L,3) + RESP_S(L,4)
      ENDDO

      RETURN
      END SUBROUTINE SOILCARB
#endif
