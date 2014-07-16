#if defined(A09_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Diagnostic RHcrit Calculation Scheme.
! Subroutine Interface:
      SUBROUTINE LS_CALC_RHCRIT(                                        &
!      Pressure related fields
     &  p_layer_centres                                                 &
!      Array dimensions
     &, LEVELS, row_length, rows, BL_LEVELS, global_row_length          &
!      Prognostic Fields
     &, T, Q, QCF, LAND, LAND_FRAC, ICE_FRAC                            &
!      Logical control
     &, l_mixing_ratio                                                  &
!      Output
     &, RHCPT)
!
      IMPLICIT NONE
!
!     Purpose: To calculate the critical relative humidity in every
!              grid-box.
!
!     Method : The critical relative humidity of a certain grid-box is
!            determined from the variance in a 3*3 set of boxes centred
!            on it. A fit, dependent on pressure, relates the variance
!            of the 3*3 region to the variance within the one grid-box.
!            This variance is converted to a critical relative humidity
!            in a straightforward fashion.
!            Some points in the 3*3 region may be excluded from the
!            variance calculation in the BL layers, if their
!            surfaces do not 'match'. The criterion for matching is that
!            land and sea-ice match, but that open ocean does not match
!            with either of these.
!            In all layers, points in the 3*3 region which lie outside 2
!            std devs of the mean are excluded in a second iteration of
!            the main calculation.
!            The best estimate of the standard deviation in a sample
!            size of n elements in which the mean is calculated from the
!            sample uses (n-1) in the denominator. This subroutine uses
!            n in the denominator because errors in the estimate of the
!            std dev from using such small sample sizes gave rise to
!            concerns over the std dev being too large occasionally,
!            and using n in the denominator rather than (n-1) addresses
!            this issue indirectly. Addressing the issue of large
!            estimates of std dev was considered more important than
!            using the 'best' estimate of the std dev.
!
! Current Owner of Code: S. Cusack / OWNER OF LS CLOUD SECTION
!
! History:
! Version   Date     Comment
!  5.1    09/03/00   Rewritten for New Dynamics.  A.C. Bushell
!                    Based on original code (part of HadAM4 package) at
!                    VN4.5 : subroutine RHCRIT_CALC.   S. Cusack
!
!  5.2    09/09/00   Changes needed to implement critical relative
!                    humidity scheme (fix polar row values).
!                    A.C. Bushell.
!  5.3    19/10/01   Use appropriate gcg routines.   S. Cusack
!
!  5.4    04/09/02   Add gcg_rvecsumf to External list. A.C. Bushell
!
!  6.1    16/08/04   New argument to subroutine added. Changes to the
!                    outlier rejection calculation.         S. Cusack
!  6.4    14/08/06   Use mixing ratio formulation.  Damian Wilson
! Code Description:
!   Language: FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
! Declarations:
!
!  Global Variables:----------------------------------------------------
!
#include "c_r_cp.h"
#include "c_epslon.h"
#include "c_lheat.h"
#include "c_0_dg_c.h"
!
!  Subroutine arguments
!-----------------------------------------------------------------------
! IN variables
!-----------------------------------------------------------------------
      INTEGER                                                           &
                        !, INTENT(IN)
     &  LEVELS                                                          &
!       No. of levels being processed.
     &, BL_LEVELS                                                       &
!       No. of BL levels being processed.
     &, row_length, rows                                                &
!       Horizontal dimensions of arrays being processed by scheme.
     &, global_row_length
!       Length of full global row (ie. counting all relevant PEs).
!
      LOGICAL                                                           &
                        !, INTENT(IN)
     &  LAND(0:row_length+1,0:rows+1)                                   &
                                          ! The model land mask
     &, l_mixing_ratio                    ! Use mixing ratios
!
      REAL                                                              &
                        !, INTENT(IN)
     &  ICE_FRAC(0:row_length+1,0:rows+1)                               &
!       Ice fraction.
     &, LAND_FRAC(0:row_length+1,0:rows+1)                              &
     &, p_layer_centres(0:row_length+1,0:rows+1,0:LEVELS)               &
!       pressure at all points, on theta levels (Pa).
!       NB: Zero level of array is surface pressure.
     &, Q(0:row_length+1,0:rows+1,LEVELS)                               &
!       Total water content, Q+QCL (kg per kg)
     &, QCF(0:row_length+1,0:rows+1,LEVELS)                             &
!       Cloud ice content at processed levels (kg water per kg air).
     &, T(0:row_length+1,0:rows+1,LEVELS)
!       Liquid/frozen water temperature (K)
!-----------------------------------------------------------------------
! OUT variables
!-----------------------------------------------------------------------
      REAL                                                              &
                        !, INTENT(OUT)
     &  RHCPT(row_length,rows,LEVELS)
!       Critical relative humidity at every gridpoint.
!
!  Local Parameters and other physical constants------------------------
!
#include "parvars.h"
#include "rhccon2b.h"
!
      REAL INV9, LS, LSRCP, ERCPR
      PARAMETER ( INV9 = 1./9.                                          &
     &          , LS = LC+LF                                            &
     &          , LSRCP = (LC+LF)/CP                                    &
     &          , ERCPR = EPSILON/(CP*R))
!
!  Local scalars--------------------------------------------------------
!
      INTEGER                                                           &
     &    I, J, K, J8                                                   &
                            ! Simple loop variables
     & ,  COUNT                                                         &
                            ! Counter of points
     & ,  IM1,IP1, JM1,JP1                                              &
                            ! (I-1),(I+1),(J-1),(J+1)
     & ,  ij_field                                                      &
                            ! Number of non-halo points in arrays
     & ,  i_length                                                      &
                            ! Row length for polar row adjustments
     & ,  i_start                                                       &
                            ! Row start point for polar row adjustments
     & ,  ISTAT             ! Status (error code) indicator
!
      REAL                                                              &
     &    MEAN_SUPSAT                                                   &
                            ! MEAN RH OF 3*3 REGION
     & ,  TOT_VAR                                                       &
                            ! TOTAL VARIANCE OF 3*3 REGION
     & ,  SUPSAT_SD_1                                                   & 
                            ! STANDARD DEVIATION OF 'R.H.' IN GRID-BOX
     & ,  SUPSAT_SD_3                                                   & 
                            ! RESOLVED STD DEV OF 'R.H.' IN 3*3 REGION
     & ,  ROOT_6                                                        &
                            ! =sqrt(6.)
     & ,  LATHT                                                         &
                            ! =Lc if T>Tm, ELSE = Lc+Lf
     & ,  TWO_SIGMA                                                     &
                          ! Two times sigma
     & ,  SUPSAT1                                                       &
                            ! Temporary variable
     & ,  SUPSAT2                                                       &
                            ! Temporary variable
     & ,  SUPSAT3                                                       &
                            ! Temporary variable
     & ,  SUPSAT4                                                       &
                            ! Temporary variable
     & ,  SUPSAT5                                                       &
                            ! Temporary variable
     & ,  SUPSAT6                                                       &
                            ! Temporary variable
     & ,  SUPSAT7                                                       &
                            ! Temporary variable
     & ,  SUPSAT8                                                       &
                            ! Temporary variable
     & ,  r_row_length      ! Reciprocal of number of points on FI row)
!
!  Local dynamic arrays-------------------------------------------------
      INTEGER                                                           &
     &   ICOUNT(row_length,rows)       ! Counter of points
!
      LOGICAL                                                           &
     &   OCEAN(0:row_length+1,0:rows+1)! Those points which are not
!                     land, and have a sea-ice fraction less than 0.25.
!
      REAL                                                              &
     &   TL(0:row_length+1,0:rows+1,LEVELS)                             &
                                            ! Conserved temperature
!                                                      (P292.1, UMDP29)
     & , QT(0:row_length+1,0:rows+1,LEVELS)                             &
                                            ! Conserved WATER
!                                                      (P292.2, UMDP29)
     & , QST(0:row_length+1,0:rows+1)                                   &
                                        ! SATURATION VAPOUR PRESSURE
     & , P_GRAD(row_length,rows,LEVELS)                                 &
                                        ! TERM WHICH RELATES
!                                         RH_SD_3 TO RH_SD_1
     & , SUPSAT(0:row_length+1,0:rows+1)                                & 
                                        ! 'RELATIVE HUMIDITY' OF GRIDBOX
     & , AL(0:row_length+1,0:rows+1)                                    &
                                        ! Defined in P292.6 in UMDP 29
     & , SURF_MULT(row_length,rows,8)                                   &
                                        ! Multiplier to take into
!                                         account surface matching.
     & , RHCPT_MEAN(LEVELS)             ! Mean of first interior row.
!
!  External subroutine calls: ------------------------------------------
      EXTERNAL QSAT
      EXTERNAL gcg_rvecsumr, gcg_rvecsumf
!- End of Header
!
! ==Main Block==--------------------------------------------------------
      ROOT_6=SQRT(6.)
!
! Levels_do1:
      DO K=1,LEVELS
! Rows_do1:
        DO J=1,rows
! Rowlen_do1:
          DO I=1,row_length
            P_GRAD(I,J,K) = p_layer_centres(I,J,K) - RHC_CON4
            P_GRAD(I,J,K) = RHC_CON1 +                                  &
     &     (RHC_CON2 * P_GRAD(I,J,K) / (RHC_CON3 + ABS(P_GRAD(I,J,K))) )
!
          END DO ! Rowlen_do1
        END DO ! Rows_do1
      END DO ! Levels_do1
!
! Ocean points defined now as not land and where ice fraction LT 0.25
! Rows_do2:
      DO J=0,rows+1
! Rowlen_do2:
        DO I=0,row_length+1
          OCEAN(I,J) = (LAND_FRAC(I,J) <  0.5) .AND.                    &
     &                                    (ICE_FRAC(I,J) <  2.5E-1)
        END DO ! Rowlen_do2
      END DO ! Rows_do2
!
! A real no. is now assigned to every neighbouring point of every point
! on the grid, if their surfaces match it has the value one, else it is
! zero.
! Eight_do1:
      DO J8=1,8
! Rows_do3:
        DO J=1,rows
! Rowlen_do3:
          DO I=1,row_length
            SURF_MULT(I,J,J8)=0.
          END DO ! Rowlen_do3
        END DO ! Rows_do3
      END DO ! Eight_do1
!
! Rows_do4:
      DO J=1,rows
        JM1 = J - 1
        JP1 = J + 1
! Rowlen_do4:
        DO I=1,row_length
          ICOUNT(I,J)=1
          IM1 = I - 1
          IP1 = I + 1
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IM1,JM1))  .OR.                  &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IM1,JM1)) ) THEN
            SURF_MULT(I,J,1) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(I,JM1))  .OR.                    &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(I,JM1)) ) THEN
            SURF_MULT(I,J,2) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IP1,JM1))  .OR.                  &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IP1,JM1)) ) THEN
            SURF_MULT(I,J,3) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IM1,J))  .OR.                    &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IM1,J)) ) THEN
            SURF_MULT(I,J,4) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IP1,J))  .OR.                    &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IP1,J)) ) THEN
            SURF_MULT(I,J,5) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IM1,JP1))  .OR.                  &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IM1,JP1)) ) THEN
            SURF_MULT(I,J,6) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(I,JP1))  .OR.                    &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(I,JP1)) ) THEN
            SURF_MULT(I,J,7) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IP1,JP1))  .OR.                  &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IP1,JP1)) ) THEN
            SURF_MULT(I,J,8) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
        END DO ! Rowlen_do4
      END DO ! Rows_do4
!
! An initial sweep is done for all grid-cells, obtaining an initial
! estimate of the variance of the 3*3 grid.
!
! Levels_do5:
      DO K=1,BL_LEVELS
! Rows_do5a:
        DO J=0,rows+1
! Rowlen_do5a:
          DO I=0,row_length+1
!   Calculate Tl and QT as in P292.1, P292.2 in UMDP 29.
!  (Assumes version 3A onwards of Section 4)
            TL(I,J,K) = T(I,J,K) - LSRCP*QCF(I,J,K)
            QT(I,J,K) = Q(I,J,K) + QCF(I,J,K)
!
          END DO ! Rowlen_do5a
        END DO ! Rows_do5a
!
! DEPENDS ON: qsat_mix
        call qsat_mix(qst(0,0),tl(0,0,K),p_layer_centres(0,0,K),        &
     &            (row_length+2)*(rows+2),l_mixing_ratio)
! Rows_do5b:
        DO J=0,rows+1
! Rowlen_do5b:
          DO I=0,row_length+1
            IF (TL(I,J,K)  >   TM) THEN
              LATHT = LC/TL(I,J,K)
            ELSE
              LATHT = LS/TL(I,J,K)
            ENDIF
            AL(I,J) = 1./ (1.+LATHT*LATHT*ERCPR*QST(I,J))
! SUPSAT given by P292.3 of UMDP 29.
            SUPSAT(I,J) = AL(I,J)*(QT(I,J,K) - QST(I,J))
          END DO ! Rowlen_do5b
        END DO ! Rows_do5b
!
! Rows_do5c:
        DO J=1,rows
          JM1 = J - 1
          JP1 = J + 1
! Rowlen_do5c:
          DO I=1,row_length
            IM1 = I - 1
            IP1 = I + 1
            SUPSAT1 = SURF_MULT(I,J,1) * SUPSAT(IM1,JM1)
            SUPSAT2 = SURF_MULT(I,J,2) * SUPSAT(I,JM1)
            SUPSAT3 = SURF_MULT(I,J,3) * SUPSAT(IP1,JM1)
            SUPSAT4 = SURF_MULT(I,J,4) * SUPSAT(IM1,J)
            SUPSAT5 = SURF_MULT(I,J,5) * SUPSAT(IP1,J)
            SUPSAT6 = SURF_MULT(I,J,6) * SUPSAT(IM1,JP1)
            SUPSAT7 = SURF_MULT(I,J,7) * SUPSAT(I,JP1)
            SUPSAT8 = SURF_MULT(I,J,8) * SUPSAT(IP1,JP1)
            MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4          &
     &                 + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8          &
     &                 + SUPSAT(I,J)) / ICOUNT(I,J)
            TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2                 &
     &              + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4                 &
     &              + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6                 &
     &              + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8                 &
     &              + SUPSAT(I,J)*SUPSAT(I,J)                           &
     &              - ICOUNT(I,J)*MEAN_SUPSAT*MEAN_SUPSAT
            TOT_VAR = ABS(TOT_VAR)
!
!  Now remove the statistical outliers from the 3*3 region, so that
!  sigma, and hence RHcrit, is not biased by extreme values.
!  Points outside 2*sigma of the mean are considered outliers and are
!  rejected.
            IF (ICOUNT(I,J)  >   1) THEN
              TWO_SIGMA = 2. * SQRT(TOT_VAR/ICOUNT(I,J))
            ELSE
              TWO_SIGMA = QST(I,J) * 0.01
            ENDIF
            COUNT=1
!
            IF (ABS(SUPSAT(IM1,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT1 = 0.
            ELSE IF (SURF_MULT(I,J,1)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(I,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT2 = 0.
            ELSE IF (SURF_MULT(I,J,2)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(IP1,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT3 = 0.
            ELSE IF (SURF_MULT(I,J,3)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(IM1,J)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT4 = 0.
            ELSE IF (SURF_MULT(I,J,4)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(IP1,J)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT5 = 0.
            ELSE IF (SURF_MULT(I,J,5)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(IM1,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT6 = 0.
            ELSE IF (SURF_MULT(I,J,6)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(I,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT7 = 0.
            ELSE IF (SURF_MULT(I,J,7)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(IP1,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT8 = 0.
            ELSE IF (SURF_MULT(I,J,8)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (COUNT >  1) THEN
              MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4        &
     &                   + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8        &
     &                   + SUPSAT(I,J)) / COUNT
              TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2               &
     &                + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4               &
     &                + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6               &
     &                + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8               &
     &                + SUPSAT(I,J)*SUPSAT(I,J)                         &
     &                - COUNT*MEAN_SUPSAT*MEAN_SUPSAT
              TOT_VAR = ABS(TOT_VAR)
              SUPSAT_SD_3 = SQRT(TOT_VAR/COUNT)
            ELSE
!           Limit the 3*3 grid variance when scatter is large.
              SUPSAT_SD_3 = QST(I,J)*0.01
            ENDIF
!
! Try to detect if the central point (i,j) is an outlier, as can happen
! in a GPS situation. If so, set the variance to a small value.
            IF (COUNT >  2) THEN
              MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4        &
     &         + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8) / (COUNT-1.0)
              TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2               &
     &                + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4               &
     &                + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6               &
     &                + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8               &
     &                - COUNT*MEAN_SUPSAT*MEAN_SUPSAT
              IF (ABS(SUPSAT(I,J)-MEAN_SUPSAT) >                        &
     &                    (2.0*SQRT(ABS(TOT_VAR)/COUNT))) THEN
                SUPSAT_SD_3 = QST(I,J)*0.01
              ENDIF
            ENDIF
!
! P_GRAD determines the relation between 3*3 and sub-grid variance.
            SUPSAT_SD_1 = P_GRAD(I,J,K) * SUPSAT_SD_3
! RHCPT defined from P292.14 in UMDP 29
            RHCPT(I,J,K) = 1. - (ROOT_6*SUPSAT_SD_1 /(AL(I,J)*QST(I,J)))
! RHcrit is now limited to lie between a range defined in RHCCON2B
            RHCPT(I,J,K) = MAX(RHCPT(I,J,K),RHC_MIN)
            RHCPT(I,J,K) = MIN(RHCPT(I,J,K),RHC_MAX)
          END DO ! Rowlen_do5c
        END DO ! Rows_do5c
!
      END DO ! Levels_do5
!
! The same calculations as above are performed, but the 'surface match'
! criterion is now dropped (atmosphere less influenced by surface at
! greater heights).
! Levels_if6:
      IF (LEVELS  >   BL_LEVELS) THEN
!
! An initial sweep is done for all grid-cells, obtaining an initial
! estimate of the variance of the 3*3 grid.
!
! Levels_do6:
        DO K=(BL_LEVELS+1),LEVELS
! Rows_do6a:
          DO J=0,rows+1
! Rowlen_do6a:
            DO I=0,row_length+1
!   Calculate Tl and QT as in P292.1, P292.2 in UMDP 29.
!  (Assumes version 3A onwards of Section 4)
              TL(I,J,K) = T(I,J,K) - LSRCP*QCF(I,J,K)
              QT(I,J,K) = Q(I,J,K) + QCF(I,J,K)
!
            END DO ! Rowlen_do6a
          END DO ! Rows_do6a
!
! DEPENDS ON: qsat_mix
          call qsat_mix(qst(0,0),tl(0,0,K),p_layer_centres(0,0,K),      &
     &              (row_length+2)*(rows+2),l_mixing_ratio)
!
! Rows_do6b:
          DO J=0,rows+1
! Rowlen_do6b:
            DO I=0,row_length+1
              IF (TL(I,J,K)  >   TM) THEN
                LATHT = LC/TL(I,J,K)
              ELSE
                LATHT = LS/TL(I,J,K)
              ENDIF
              AL(I,J) = 1./ (1.+LATHT*LATHT*ERCPR*QST(I,J))
! SUPSAT given by P292.3 of UMDP 29.
              SUPSAT(I,J) = AL(I,J)*(QT(I,J,K) - QST(I,J))
            END DO ! Rowlen_do6b
          END DO ! Rows_do6b
!
! Rows_do6c:
          DO J=1,rows
            JM1 = J - 1
            JP1 = J + 1
! Rowlen_do6c:
            DO I=1,row_length
              IM1 = I - 1
              IP1 = I + 1
              SUPSAT1 = SUPSAT(IM1,JM1)
              SUPSAT2 = SUPSAT(I,JM1)
              SUPSAT3 = SUPSAT(IP1,JM1)
              SUPSAT4 = SUPSAT(IM1,J)
              SUPSAT5 = SUPSAT(IP1,J)
              SUPSAT6 = SUPSAT(IM1,JP1)
              SUPSAT7 = SUPSAT(I,JP1)
              SUPSAT8 = SUPSAT(IP1,JP1)
              MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4        &
     &                   + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8        &
     &                   + SUPSAT(I,J)) * INV9
              TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2               &
     &                + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4               &
     &                + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6               &
     &                + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8               &
     &                + SUPSAT(I,J)*SUPSAT(I,J)                         &
     &                - 9. * MEAN_SUPSAT*MEAN_SUPSAT
              TOT_VAR = ABS(TOT_VAR)
!
!  Now remove the statistical outliers from the 3*3 region, so that
!  sigma, and hence RHcrit, is not biased by extreme values.
!  Points outside 2*sigma of the mean are considered outliers and are
!  rejected.
              TWO_SIGMA = 0.67*SQRT(TOT_VAR)  ! =2*SQRT(TOT_VAR/9)
              COUNT=1
!
            IF (ABS(SUPSAT(IM1,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT1=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(I,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT2=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(IP1,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT3=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(IM1,J)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT4=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(IP1,J)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT5=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(IM1,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT6=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(I,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT7=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(IP1,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT8=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
              IF (COUNT  >   1) THEN
                MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4      &
     &                     + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8      &
     &                     + SUPSAT(I,J)) / COUNT
                TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2             &
     &                  + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4             &
     &                  + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6             &
     &                  + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8             &
     &                  + SUPSAT(I,J)*SUPSAT(I,J)                       &
     &                  - COUNT*MEAN_SUPSAT*MEAN_SUPSAT
                TOT_VAR = ABS(TOT_VAR)
                SUPSAT_SD_3 = SQRT(TOT_VAR/COUNT)
              ELSE
!             Limit the 3*3 grid variance when scatter is large.
                SUPSAT_SD_3 = QST(I,J) * 0.01
              ENDIF
!
! Try to detect if the central point (i,j) is an outlier, as can happen
! in a GPS situation. If so, set the variance to a small value.
            IF (COUNT >  2) THEN
              MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4        &
     &         + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8) / (COUNT-1.0)
              TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2               &
     &                + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4               &
     &                + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6               &
     &                + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8               &
     &                - COUNT*MEAN_SUPSAT*MEAN_SUPSAT
              IF (ABS(SUPSAT(I,J)-MEAN_SUPSAT) >                        &
     &                    (2.0*SQRT(ABS(TOT_VAR)/COUNT))) THEN
                SUPSAT_SD_3 = QST(I,J)*0.01
              ENDIF
            ENDIF
!
! P_GRAD determines the relation between 3*3 and sub-grid variance.
              SUPSAT_SD_1 = P_GRAD(I,J,K) * SUPSAT_SD_3
! RHCPT defined from P292.14 in UMDP 29
              RHCPT(I,J,K) = 1.-(ROOT_6*SUPSAT_SD_1 /(AL(I,J)*QST(I,J)))
! RHcrit is now limited to lie between a range defined in RHCCON2B
              RHCPT(I,J,K) = MAX(RHCPT(I,J,K),RHC_MIN)
              RHCPT(I,J,K) = MIN(RHCPT(I,J,K),RHC_MAX)
            END DO ! Rowlen_do6c
          END DO ! Rows_do6c
!
        END DO  ! Levels_do6
      ENDIF  ! Levels_if6
!
! Tidy up at South Pole : Pole is mean of first interior row.
! SouthPole_if1:
      IF (at_extremity(PSouth)) THEN
!
        ij_field = row_length * rows
        i_length = row_length
! Start point of first interior (ie. non-polar) row
        i_start  = row_length + 1
!
        r_row_length = 1. / global_row_length
!       Number of points in sum should be global_row_length: might
!       need to modify r_row_length if i_length lt row_length.
!       r_row_length = r_row_length * row_length / i_length
!
! Sum over points in PEs in order along first interior row
! (gc_proc_row_group is group ID for rows of PEs, here only PSouth).
#if defined(REPROD)
        CALL gcg_rvecsumr(ij_field, i_length, i_start, LEVELS, RHCPT,   &
     &                    gc_proc_row_group, ISTAT, RHCPT_MEAN)
#else
        CALL gcg_rvecsumf(ij_field, i_length, i_start, LEVELS, RHCPT,   &
     &                    gc_proc_row_group, ISTAT, RHCPT_MEAN)
#endif
!
! Levels_do7:
        DO K=1,LEVELS
          RHCPT_MEAN(K) = RHCPT_MEAN(K) * r_row_length
! Rowlen_do7:
          DO I=1,row_length
            RHCPT(I,1,K) = RHCPT_MEAN(K)
          END DO ! Rowlen_do7
        END DO  ! Levels_do7
!
      ENDIF  ! SouthPole_if1
!
! Tidy up at North Pole : Pole is mean of first interior row.
! NorthPole_if1:
      IF (at_extremity(PNorth)) THEN
!
        ij_field = row_length * rows
        i_length = row_length
! Start point of first interior (ie. non-polar) row
        i_start  = (rows - 2) * row_length + 1
!
        r_row_length = 1. / global_row_length
!       Number of points in sum should be global_row_length: might
!       need to modify r_row_length if i_length lt row_length.
!       r_row_length = r_row_length * row_length / i_length
!
! Sum over points in PEs in order along first interior row
! (gc_proc_row_group is group ID for rows of PEs, here only PNorth).
#if defined(REPROD)
        CALL gcg_rvecsumr(ij_field, i_length, i_start, LEVELS, RHCPT,   &
     &                    gc_proc_row_group, ISTAT, RHCPT_MEAN)
#else
        CALL gcg_rvecsumf(ij_field, i_length, i_start, LEVELS, RHCPT,   &
     &                    gc_proc_row_group, ISTAT, RHCPT_MEAN)
#endif
!
! Levels_do8:
        DO K=1,LEVELS
          RHCPT_MEAN(K) = RHCPT_MEAN(K) * r_row_length
! Rowlen_do8:
          DO I=1,row_length
            RHCPT(I,rows,K) = RHCPT_MEAN(K)
          END DO ! Rowlen_do8
        END DO  ! Levels_do8
!
      ENDIF  ! NorthPole_if1
!
      RETURN
      END SUBROUTINE LS_CALC_RHCRIT
#endif
