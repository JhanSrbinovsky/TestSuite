#if defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine COMPETE ------------------------------------------------
!!!
!!! Purpose : Updates fractional coverage of each functional type.
!!!           Requires a dominance hierachy as input.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  4.4     10/97     New Deck. Peter Cox
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.
!!!                    Richard Betts
!!!  5.2   15/11/00    Re-written for New Dynamics     M. Best
!!!
!!!END ----------------------------------------------------------------
      SUBROUTINE COMPETE (DOM,LAND_PTS,TRIF_PTS,TRIF_INDEX              &
     &,                   B,DB_DFRAC,FORW,GAMMA,NOSOIL                  &
     &,                   FRAC,DFRAC)

      IMPLICIT NONE

#include "nstypes.h"

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,TRIF_PTS                                                         &
                                  ! IN Number of points on which
!                                 !    TRIFFID may operate
     &,K,L,M,N,T                  ! WORK Loop counters.

      INTEGER                                                           &
     & DOM(LAND_PTS,NPFT)                                               &
                                  ! IN Dominance hierachy.
     &,TRIF_INDEX(LAND_PTS)       ! IN Indices of land points on
!                                 !    which TRIFFID may operate

      REAL                                                              &
     & B(LAND_PTS,NPFT)                                                 &
                                  ! IN Mean rate of change of
!                                 !    vegetation fraction over
!                                 !    the timestep (kg C/m2/360days).
     &,DB_DFRAC(LAND_PTS,NPFT,NPFT)                                     &
!                                 ! IN Rate of change of B
!                                 !    with vegetation fraction.
     &,FORW                                                             &
                                  ! IN Forward weighting factor.
     &,GAMMA                                                            &
                                  ! IN Inverse timestep (/360days).
     &,NOSOIL(LAND_PTS)                                                 &
                                  ! IN Fractional area not available
!                                 !    to vegetation.
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                  ! INOUT Updated areal fraction.
     &,DFRAC(LAND_PTS,NPFT)                                             &
                                  ! OUT Increment to areal fraction.
     &,DENOM                                                            &
                                  ! WORK Denominator of update
!                                 !      equation.
     &,FRACN,FRACM                                                      &
                                  ! WORK Fractions used in the spreading
!                                 !      calculation.
     &,NUMER                                                            &
                                  ! WORK Numerator of the update
!                                 !      equation.
     &,SPACE(LAND_PTS)                                                  &
                                  ! WORK Available space.
     &,P1,P2,Q1,Q2,R1,R2          ! WORK Coefficients in simultaneous
!                                 !      equations.
!----------------------------------------------------------------------
! Local parameters
!----------------------------------------------------------------------
#include "seed.h"
#include "descent.h"


!----------------------------------------------------------------------
! Initialisations. Set increments to zero and define the space
! available to the dominant type leaving space for the seeds of others.
!----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)
        DO N=1,NPFT
          DFRAC(L,N) = 0.0
        ENDDO
        SPACE(L) = 1-NOSOIL(L)-FRAC_MIN*(NPFT-1)
      ENDDO

!----------------------------------------------------------------------
! Calculate the increments to the tree fractions
!----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)
        N = DOM(L,1)
        M = DOM(L,2)

        FRACN=FRAC(L,N)
        FRACN=MAX(FRACN,FRAC_SEED)

        FRACM=FRAC(L,M)
        FRACM=MAX(FRACM,FRAC_SEED)

        P1 = GAMMA/FRACN-FORW*DB_DFRAC(L,N,N)
        P2 = GAMMA/FRACM-FORW*DB_DFRAC(L,M,M)
        Q1 = -FORW*DB_DFRAC(L,N,M)
        Q2 = -FORW*DB_DFRAC(L,M,N)
        R1 = B(L,N)
        R2 = B(L,M)
        DO K=1,NPFT
          R1 = R1+FORW*(DB_DFRAC(L,N,K)*DFRAC(L,K))
          R2 = R2+FORW*(DB_DFRAC(L,M,K)*DFRAC(L,K))
        ENDDO

        NUMER = R1-(Q1/P2)*R2
        DENOM = P1-(Q1/P2)*Q2
        DENOM = MAX(DENOM,DENOM_MIN)
        DFRAC(L,N) = NUMER/DENOM
        FRAC(L,N) = FRAC(L,N)+DFRAC(L,N)

        IF (FRAC(L,N) <  FRAC_MIN) THEN
          DFRAC(L,N) = DFRAC(L,N)+(FRAC_MIN-FRAC(L,N))
          FRAC(L,N) = FRAC_MIN
        ELSEIF (FRAC(L,N) >  SPACE(L)) THEN
          DFRAC(L,N) = DFRAC(L,N)+(SPACE(L)-FRAC(L,N))
          FRAC(L,N) = SPACE(L)
        ENDIF

        SPACE(L) = SPACE(L)-FRAC(L,N)+FRAC_MIN

        NUMER = R2-Q2*DFRAC(L,N)
        DENOM = P2
        DENOM = MAX(DENOM,DENOM_MIN)
        DFRAC(L,M) = NUMER/DENOM
        FRAC(L,M) = FRAC(L,M)+DFRAC(L,M)

        IF (FRAC(L,M) <  FRAC_MIN) THEN
          DFRAC(L,M) = DFRAC(L,M)+(FRAC_MIN-FRAC(L,M))
          FRAC(L,M) = FRAC_MIN
        ELSEIF (FRAC(L,M) >  SPACE(L)) THEN
          DFRAC(L,M) = DFRAC(L,M)+(SPACE(L)-FRAC(L,M))
          FRAC(L,M) = SPACE(L)
        ENDIF

        SPACE(L) = SPACE(L)-FRAC(L,M)+FRAC_MIN

      ENDDO

!----------------------------------------------------------------------
! Calculate the increment to the shrub fraction
!----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)
        N = DOM(L,3)

        FRACN=FRAC(L,N)
        FRACN=MAX(FRACN,FRAC_SEED)

        DENOM = GAMMA/FRACN-FORW*DB_DFRAC(L,N,N)
        DENOM = MAX(DENOM,DENOM_MIN)

        NUMER = B(L,N)
        DO K=1,NPFT
          NUMER = NUMER+FORW*(DB_DFRAC(L,N,K)*DFRAC(L,K))
        ENDDO

        DFRAC(L,N) = NUMER/DENOM
        FRAC(L,N) = FRAC(L,N)+DFRAC(L,N)

        IF (FRAC(L,N) <  FRAC_MIN) THEN
          DFRAC(L,N) = DFRAC(L,N)+(FRAC_MIN-FRAC(L,N))
          FRAC(L,N) = FRAC_MIN
        ELSEIF (FRAC(L,N) >  SPACE(L)) THEN
          DFRAC(L,N) = DFRAC(L,N)+(SPACE(L)-FRAC(L,N))
          FRAC(L,N) = SPACE(L)
        ENDIF

        SPACE(L) = SPACE(L)-FRAC(L,N)+FRAC_MIN
      ENDDO


!----------------------------------------------------------------------
! Calculate the increments to the grass fractions
!----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)
        N = DOM(L,4)
        M = DOM(L,5)

        FRACN=FRAC(L,N)
        FRACN=MAX(FRACN,FRAC_SEED)

        FRACM=FRAC(L,M)
        FRACM=MAX(FRACM,FRAC_SEED)

        P1 = GAMMA/FRACN-FORW*DB_DFRAC(L,N,N)
        P2 = GAMMA/FRACM-FORW*DB_DFRAC(L,M,M)
        Q1 = -FORW*DB_DFRAC(L,N,M)
        Q2 = -FORW*DB_DFRAC(L,M,N)
        R1 = B(L,N)
        R2 = B(L,M)
        DO K=1,NPFT
          R1 = R1+FORW*(DB_DFRAC(L,N,K)*DFRAC(L,K))
          R2 = R2+FORW*(DB_DFRAC(L,M,K)*DFRAC(L,K))
        ENDDO

        NUMER = R1-(Q1/P2)*R2
        DENOM = P1-(Q1/P2)*Q2
        DENOM = MAX(DENOM,DENOM_MIN)
        DFRAC(L,N) = NUMER/DENOM
        FRAC(L,N) = FRAC(L,N)+DFRAC(L,N)

        IF (FRAC(L,N) <  FRAC_MIN) THEN
          DFRAC(L,N) = DFRAC(L,N)+(FRAC_MIN-FRAC(L,N))
          FRAC(L,N) = FRAC_MIN
        ELSEIF (FRAC(L,N) >  SPACE(L)) THEN
          DFRAC(L,N) = DFRAC(L,N)+(SPACE(L)-FRAC(L,N))
          FRAC(L,N) = SPACE(L)
        ENDIF

        SPACE(L) = SPACE(L)-FRAC(L,N)+FRAC_MIN

        NUMER = R2-Q2*DFRAC(L,N)
        DENOM = P2
        DENOM = MAX(DENOM,DENOM_MIN)
        DFRAC(L,M) = NUMER/DENOM
        FRAC(L,M) = FRAC(L,M)+DFRAC(L,M)

        IF (FRAC(L,M) <  FRAC_MIN) THEN
          DFRAC(L,M) = DFRAC(L,M)+(FRAC_MIN-FRAC(L,M))
          FRAC(L,M) = FRAC_MIN
        ELSEIF (FRAC(L,M) >  SPACE(L)) THEN
          DFRAC(L,M) = DFRAC(L,M)+(SPACE(L)-FRAC(L,M))
          FRAC(L,M) = SPACE(L)
        ENDIF

        SPACE(L) = SPACE(L)-FRAC(L,M)+FRAC_MIN

      ENDDO

!----------------------------------------------------------------------
! Diagnose the new bare soil fraction
!----------------------------------------------------------------------
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)
        FRAC(L,SOIL) = 1.0-NOSOIL(L)
        DO N=1,NPFT
          FRAC(L,SOIL) = FRAC(L,SOIL)-FRAC(L,N)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE COMPETE
#endif
