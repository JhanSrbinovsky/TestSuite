#if defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine DECAY --------------------------------------------------
!!!
!!! Purpose : Updates carbon contents of the soil.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4     10/97     New Deck. Peter Cox
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.
!!!                    Richard Betts
!!!  5.2   15/11/00    Re-written for New Dynamics      M. Best
!  6.2  01/03/06  extend to incorporate RothC.              C.D. Jones
!!!
!!!END ----------------------------------------------------------------
      SUBROUTINE DECAY (LAND_PTS,TRIF_PTS,TRIF_INDEX                    &
     &,                 DPC_DCS,FORW,GAMMA,PC,CS)

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
     &,L,T,N                      ! WORK Loop counters

      REAL                                                              &
     & DPC_DCS(LAND_PTS,4)                                              &
                                  ! IN Rate of change of PC with
!                                 !    soil carbon (yr).
     &,FORW                                                             &
                                  ! IN Forward timestep weighting.
     &,GAMMA                                                            &
                                  ! IN Inverse timestep (/360days).
     &,PC(LAND_PTS,4)                                                   &
                                  ! IN Net carbon flux into the
!                                 !    soil (kg C/m2/360days).
     &,CS(LAND_PTS,4)                                                   &
                                  ! INOUT Soil carbon (kg C/m2).
     &,DENOM                                                            &
                                  ! WORK Denominator of update
!                                 !      equation.
     &,NUMER                      ! WORK Numerator of the update
!                                 !      equation.
!----------------------------------------------------------------------
! Local parameters
!----------------------------------------------------------------------
#include "csmin.h"
#include "descent.h"

      DO N=1,4
        DO T=1,TRIF_PTS
          L=TRIF_INDEX(T)

          NUMER = PC(L,N)
          DENOM = GAMMA+FORW*DPC_DCS(L,N)
          DENOM = MAX(DENOM,DENOM_MIN)

          CS(L,N) = CS(L,N)+NUMER/DENOM

          CS(L,N) = MAX(CS_MIN,CS(L,N))

        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DECAY
#endif
