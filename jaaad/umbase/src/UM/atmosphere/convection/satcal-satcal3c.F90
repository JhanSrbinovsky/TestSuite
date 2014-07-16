#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SATCAL-------------------------------------------------
!LL
!LL  PURPOSE : CALCULATES SATURATED TEMPERATURE
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL   4.4   17/10/97  New version optimised for T3E.
!LL                   Single PE optimisations           D.Salmond
!LL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!     6.2   03/02/05  Added section 5A. R A Stratton
!     6.4   12/12/06  Removed def 3C. R A Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL  VERSION NO. 4  DATED 5/2/92
!LL
!LL  SYSTEM TASK : P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE SATCAL (NPNTS,T,TH,PK,QS,THDDS,K,EXK,Q_K,THE_K)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS
!-----------------------------------------------------------------------
!
#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_0_dg_c.h"
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
!
      INTEGER I                 ! LOOP COUNTER
!
      INTEGER IC                ! LOOP COUNTER
!
      INTEGER NPNTS             ! VECTOR LENGTH
!
      INTEGER K                 ! IN PRESENT MODEL LAYER
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
      REAL TH(NPNTS)            ! IN POTENTIAL TEMPERATURE (K)
!
      REAL T(NPNTS)             ! IN TEMPERATURE (K)
!
      REAL PK(NPNTS)            ! IN PRESSURE OF LAYER K (PA)
!
      REAL Q_K(NPNTS)           ! IN MIXING RATIO OF LAYER K (KG/KG)
!
      REAL EXK(NPNTS)           ! IN EXNER RATIO OF LAYER K
!
      REAL THE_K(NPNTS)         ! IN ENVIRONMENTAL POTENTIAL TEMPERATURE
                                !    IN LAYER K
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!-----------------------------------------------------------------------
!
      REAL QS(NPNTS)            ! OUT SATURATED SPECIFIC HUMIDITY
                                !     (KG/KG)
!
      REAL THDDS(NPNTS)         ! OUT SATURATED ENVIRONMENTAL
                                !     POTENTIAL TEMPERATURE (K)
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE LOCALLY DEFINED
!-----------------------------------------------------------------------
!
      REAL TS(NPNTS)            ! SATURATED TEMPERATURE (K)
!
      REAL T_FG(NPNTS)          ! TEMPERATURE FIRST GUESS (K)
!
      REAL TH_FG(NPNTS)         ! POTENTIAL TEMPERATURE FIRST GUESS (K)
!
      REAL DQBYDT(NPNTS)        ! FIRST GUESS AT MIXING RATIO INCREMENT
                                ! (KG/KG/K)
!
      REAL L                    ! LATENT HEAT
!
!
!-----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!-----------------------------------------------------------------------
!
      EXTERNAL QSAT, DQS_DTH
!
!-----------------------------------------------------------------------
! SET INITIAL FIRST GUESS TEMPERATURE AND THETA - BASED UPON
! ENVIRONMENTAL TEMPERATURE IN LAYER K
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
       TH_FG(I) = THE_K(I)
       T_FG(I) = TH_FG(I)*EXK(I)
      END DO
!
!----------------------------------------------------------------------
! CALCULATE QSAT FOR INITIAL FIRST GUESS TEMPERATURE
!----------------------------------------------------------------------
!
! DEPENDS ON: qsat
      CALL QSAT(QS,T_FG,PK,NPNTS)
!
!----------------------------------------------------------------------
! DO TWO ITERATIONS TO FIND SATURATION POINT DUE TO EVAPORATION
!----------------------------------------------------------------------
!
      DO IC=1,2
!
!----------------------------------------------------------------------
! CALCULATE DQSAT/DT FOR FIRST GUESS TEMPERATURE
!----------------------------------------------------------------------
!
! DEPENDS ON: dqs_dth
       CALL DQS_DTH(DQBYDT,K,TH_FG,QS,EXK,NPNTS)
!
!----------------------------------------------------------------------
! CALCULATE UPDATED TEMPERATURE AT SATURATION
!----------------------------------------------------------------------
!
       DO I=1,NPNTS
!
        IF (T_FG(I) >  TM) THEN
         L=LC
        ELSE
         L=LC+LF
        END IF
!
        THDDS(I) = (TH(I)*CP*EXK(I) - L*(QS(I)-Q_K(I)-                  &
     &                  TH_FG(I)*DQBYDT(I))) /                          &
     &                  (CP*EXK(I)+L*DQBYDT(I))
!
!----------------------------------------------------------------------
! CALCULATE TEMPERATURE AT SATURATION AND UPDATE FIRST GUESS
!----------------------------------------------------------------------
!
        TH_FG(I) = THDDS(I)
        T_FG(I) = TH_FG(I)*EXK(I)
!
       END DO
!
!----------------------------------------------------------------------
! CALCULATE REVISED SATURATION MIXING RATIO AT SATURATION
!---------------------------------------------------------------------
!
! DEPENDS ON: qsat
       CALL QSAT(QS,T_FG,PK,NPNTS)
!
      END DO
!
      RETURN
      END SUBROUTINE SATCAL
!
#endif
