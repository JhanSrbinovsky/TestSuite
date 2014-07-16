#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE BOUY_TQ

! PURPOSE: To calculate buoyancy parameters on p,T,q-levels
!
! METHOD:
!
! HISTORY:
! DATE   VERSION   COMMENT
! ----   -------   -------
! new deck
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!LL   5.2  27/09/00   change from QSAT_2D to QSAT           A.Malcolm
!!!  5.3    Feb. 01  Calculate grid-box mean parameters (APLock)
!  6.1  17/05/04  Change Q, QCL, QCF dims to enable substepping.
!                                                       M. Diamantakis
!
! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
!   THIS CODE IS WRITTEN TO UMDP 3 PROGRAMMING STANDARDS.
!

      SUBROUTINE BOUY_TQ (                                              &
     & row_length, rows, halo_i,halo_j,BL_LEVELS, LQ_MIX_BL             &
     &,P,T,Q,QCF,QCL,CF                                                 &
     &,BT,BQ,BT_CLD,BQ_CLD,BT_GB,BQ_GB,A_QS,A_DQSDT,DQSDT               &
     &,LTIMER                                                           &
     & )

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL LTIMER          ! IN Flag for TIMER diagnostics


      LOGICAL LQ_MIX_BL       ! IN switch for using mixing ratios
      INTEGER                                                           &
     & row_length, rows, halo_i,halo_j                                  &
     &,BL_LEVELS              ! IN No. of atmospheric levels for which
!                                boundary layer fluxes are calculated.
!                                Assumed  <=30 for dimensioning GAMMA()
!                                in common deck C_GAMMA

      REAL                                                              &
     & P(row_length,rows,BL_LEVELS)                                     &
                                      ! IN Pressure at pressure points.
     &,T(row_length,rows,BL_LEVELS)                                     &
                                      ! IN Temperature (K). At P points
     &,Q(row_length,rows,BL_LEVELS)                                     &
                                ! IN Sp humidity (kg water per kg air).
     &,QCL(row_length,rows,BL_LEVELS)                                   &
                                ! IN Cloud liq water (kg per kg air).
     &,QCF(row_length,rows,BL_LEVELS)                                   &
                                ! IN Cloud liq water (kg per kg air).
     &,CF(row_length, rows, BL_LEVELS)! IN Cloud fraction (decimal).


! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

      REAL                                                              &
     & BQ(row_length,rows,BL_LEVELS)                                    &
                                ! OUT A buoyancy parameter for clear air
     &,BT(row_length,rows,BL_LEVELS)                                    &
                                ! OUT A buoyancy parameter for clear air
     &,BQ_CLD(row_length,rows,BL_LEVELS)                                &
!                             ! OUT A buoyancy parameter for cloudy air
     &,BT_CLD(row_length,rows,BL_LEVELS)                                &
!                             ! OUT A buoyancy parameter for cloudy air
     &,BQ_GB(row_length,rows,BL_LEVELS)                                 &
                                ! OUT A grid-box mean buoyancy parameter
     &,BT_GB(row_length,rows,BL_LEVELS)                                 &
                                ! OUT A grid-box mean buoyancy parameter
     &,A_QS(row_length,rows,BL_LEVELS)                                  &
!                             ! OUT Saturated lapse rate factor
     &,A_DQSDT(row_length,rows,BL_LEVELS)                               &
!                             ! OUT Saturated lapse rate factor
     &,DQSDT(row_length,rows,BL_LEVELS)
!                             ! OUT Derivative of q_SAT w.r.t. T

! LOCAL VARIABLES.

      REAL                                                              &
     & QS(row_length,rows)            ! WORK Saturated mixing ratio.

      INTEGER                                                           &
     &  I,j                                                             &
     &, K

      REAL                                                              &
     &  BC

      EXTERNAL                                                          &
     &  QSAT, TIMER

#include "c_0_dg_c.h"
#include "c_lheat.h"
#include "c_g.h"
#include "c_r_cp.h"
#include "c_epslon.h"
#include "c_vkman.h"
#include "c_soilh.h"
#include "c_mdi.h"


      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP
      PARAMETER (                                                       &
     & ETAR=1.0/(1.0-EPSILON)                                           &
                                ! Used in buoyancy parameter BETAC.
     &,GRCP=G/CP                                                        &
                                ! Used in DZTL, FTL calculations.
     &,LCRCP=LC/CP                                                      &
                                ! Latent heat of condensation / CP.
     &,LFRCP=LF/CP                                                      &
                                ! Latent heat of fusion / CP.
     &,LS=LC+LF                                                         &
                                ! Latent heat of sublimation.
     &,LSRCP=LS/CP                                                      &
                                ! Latent heat of sublimation / CP.
     &)

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BOUY_TQ ',3)
      ENDIF
!-----------------------------------------------------------------------
!! 1.  Loop round levels.
!-----------------------------------------------------------------------
      DO K=1,BL_LEVELS
!-----------------------------------------------------------------------
!! 1.1 Calculate saturated specific humidity at pressure and
!!     temperature of current level.
!-----------------------------------------------------------------------
! DEPENDS ON: qsat_mix
        CALL QSAT_mix(QS,T(1,1,K),P(1,1,K),row_length*rows,Lq_mix_bl)
!
        do j=1,rows
        DO I=1,row_length
!ajm        DO I=P1,P1+P_POINTS-1

!-----------------------------------------------------------------------
!! 1.2 Calculate buoyancy parameters BT and BQ, required for the
!!     calculation of stability.
!-----------------------------------------------------------------------

          BT(I,j,K) = 1.0/T(I,j,K)
          BQ(I,j,K) =                                                   &
     &      C_VIRTUAL/(1.0+C_VIRTUAL*Q(I,j,K)-QCL(I,j,K)-QCF(I,j,K))
!
          IF (T(I,j,K)  >   TM) THEN
            DQSDT(I,j,K) = (EPSILON * LC * QS(I,j))                     &
     &                   / ( R * T(I,j,K) * T(I,j,K) )
!                      ...  (Clausius-Clapeyron) for T above freezing
!
            A_QS(I,j,K) = 1.0 / (1.0 + LCRCP*DQSDT(I,j,K))
!
            A_DQSDT(I,j,K) = A_QS(I,j,K) * DQSDT(I,j,K)
!
            BC = LCRCP*BT(I,j,K) - ETAR*BQ(I,j,K)
!
          ELSE
            DQSDT(I,j,K) = (EPSILON * LS * QS(I,j))                     &
     &                   / ( R * T(I,j,K) * T(I,j,K) )
!                      ...  (Clausius-Clapeyron) for T below freezing
!
            A_QS(I,j,K) = 1.0 / (1.0 + LSRCP*DQSDT(I,j,K))
!
            A_DQSDT(I,j,K) = A_QS(I,j,K) * DQSDT(I,j,K)
!
            BC = LSRCP*BT(I,j,K) - ETAR*BQ(I,j,K)
!
          ENDIF
!
!-----------------------------------------------------------------------
!! 1.3 Calculate in-cloud buoyancy parameters.
!-----------------------------------------------------------------------
!
          BT_CLD(I,j,K) = BT(I,j,K) - A_DQSDT(I,j,K) * BC
          BQ_CLD(I,j,K) = BQ(I,j,K) + A_QS(I,j,K) * BC

!-----------------------------------------------------------------------
!! 1.4 Calculate grid-box mean buoyancy parameters.
!-----------------------------------------------------------------------
!
          BT_GB(I,j,K) = BT(I,j,K) +                                    &
     &                   CF(I,j,K)*( BT_CLD(I,j,K) - BT(I,j,K) )
          BQ_GB(I,j,K) = BQ(I,j,K) +                                    &
     &                   CF(I,j,K)*( BQ_CLD(I,j,K) - BQ(I,j,K) )
!
        ENDDO ! p_points,j
        ENDDO ! p_points,i
      ENDDO ! bl_levels

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BOUY_TQ ',4)
      ENDIF
      RETURN
      END SUBROUTINE BOUY_TQ
#endif
