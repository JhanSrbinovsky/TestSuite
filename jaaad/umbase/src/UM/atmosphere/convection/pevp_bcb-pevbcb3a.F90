#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PEVP_BCB-----------------------------------------------
!LL
!LL  PURPOSE : EVAPORATE RAIN BELOW CLOUD BASE IF NO DOWNDRAUGHT
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL   4.0   5/05/95   New deck at version 4.0 to include a pressure
!LL                   dependency in the calculation of evaporation
!LL                   of convective precipitation.
!LL                   Pete Inness.
!LL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!     5.5  17/04/03   Removal of reference to obsolete section
!                     A05_3B. T.White
!     6.2  03/02/05  Added section 5A. R A Stratton
!     6.4  12/12/06  Removed def 3C. R A Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL  VERSION NO. 4 DATED 5/2/92
!LL
!LL  LOGICAL COMPONENTS COVERED:
!LL
!LL  SYSTEM TASK : P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE PEVP_BCB (NPNTS,K,ICCB,TH,PK,Q,DELP,RAIN,SNOW,         &
     &                     DTHBYDT,DQBYDT,EXK,TIMESTEP,CCA)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!  CONSTANTS
!-----------------------------------------------------------------------
!
#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_g.h"
#include "cldarea.h"
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
!
      INTEGER I                  ! IN LOOP COUNTER
!
      INTEGER NPNTS              ! VECTOR LENGTH
!
      INTEGER K                  ! IN PRESENT MODEL LAYER
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
      INTEGER ICCB(NPNTS)        ! IN CONVECTIVE CLOUD BASE LAYER
!
      REAL PK(NPNTS)             ! IN PRESSURE (PA)
!
      REAL Q(NPNTS)              ! IN MIXING RATIO (KG/KG)
!
      REAL TH(NPNTS)             ! IN POTENTIAL TEMPERATURE (K)
!
      REAL DELP(NPNTS)           ! IN CHANGE IN PRESSURE ACROSS
                                 !    LAYER K-1 (PA)
!
      REAL EXK(NPNTS)            ! IN EXNER RATIO OF LAYER K
!
      REAL TIMESTEP              ! IN MODEL TIMESTEP (S)
!
      REAL CCA(NPNTS)            ! IN CONVECTIVE CLOUD AMOUNT
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!-----------------------------------------------------------------------
!
      REAL DTHBYDT(NPNTS)        ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE (K/S)
!
      REAL DQBYDT(NPNTS)         ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING RATIO
                                 !     (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MIXING RATIO
                                 !     AFTER EVAPORATION BELOW CLOUD
                                 !     BASE (KG/KG/S)
!
      REAL RAIN(NPNTS)           ! INOUT
                                 ! IN  AMOUNT OF FALLING RAIN
                                 !     (KG/M**2/S)
                                 ! OUT UPDATED AMOUNT OF FALLING RAIN
                                 !     (KG/M**2/S)
!
      REAL SNOW(NPNTS)           ! INOUT
                                 ! IN  AMOUNT OF FALLING SNOW
                                 !     (KG/M**2/S)
                                 ! OUT UPDATED AMOUNT OF FALLING SNOW
                                 !     (KG/M**2/S)
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
!
      REAL T(NPNTS)              ! MODEL TEMPERATURE (K)
!
      REAL EVAP_RAIN(NPNTS)      ! AMOUNT OF EVAPORATION OF RAIN
!
      REAL SUB_SNOW(NPNTS)       ! AMOUNT OF SNOW SUBLIMATION
!
      REAL QSATE(NPNTS)          ! SATURATED MIXING RATIO IN
                                 ! ENVIRONMENT (KG/KG)
!
      REAL DELQ(NPNTS)           ! CHANGE IN MIXING RATIO ACROSS LAYER K
                                 ! (KG/KG)
!
      REAL THS(NPNTS)            ! SATURATED PARCEL POTENTIAL
                                 ! TEMPERATURE (K)
!
      REAL QS(NPNTS)             ! SATURATED PARCEL MIXING RATIO
!
      LOGICAL BEVAP(NPNTS)       ! MASK FOR THOSE POINTS WHERE
                                 ! EVAPORATION OCCURS
!
      REAL DTHBYDT_EVP(NPNTS)    ! INCREMENT TO POTENTIAL TEMPERATURE
                                 ! DUE TO EVAPORATION (K)
!
      REAL DQBYDT_EVP(NPNTS)     ! INCREMENT TO MIXING RATIO DUE TO
                                 ! EVAPORATION (KG/KG)
!
      REAL DTHBYDT_SAT(NPNTS)    ! INCREMENT TO POTENTIAL TEMPERATURE
                                 ! DUE TO SATURATION (K)
!
      REAL FACTOR(NPNTS)         ! DTHBYDT_SAT / DTHBYDT_EVP
!
      REAL RHO(NPNTS)            ! DENSITY OF AIR IN PARCEL
!
!-----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!-----------------------------------------------------------------------
!
      EXTERNAL QSAT, EVP, SATCAL
!
!-----------------------------------------------------------------------
! EVAPORATE RAIN IN LAYER K IF LAYER K IS BELOW CLOUD BASE
! CALCULATE MOISTURE SUB-SATURATION
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
        T(I) = TH(I)*EXK(I)
        BEVAP(I) = .FALSE.
      END DO
!
! DEPENDS ON: qsat
      CALL QSAT(QSATE,T,PK,NPNTS)
!
      DO I=1,NPNTS
       IF (K  <   ICCB(I)) THEN
         DELQ(I) = QSATE(I)-Q(I)
!
!-----------------------------------------------------------------------
! CHECK IF EVAPORATION POSSIBLE
!-----------------------------------------------------------------------
!
         IF ((RAIN(I) >  0.0 .OR. SNOW(I) >  0.0) .AND.                 &
     &        DELQ(I)  >   0.0) THEN
!
            BEVAP(I) = .TRUE.
            RHO(I) = PK(I) / (R*T(I))
         END IF
       END IF
      END DO
!
!-----------------------------------------------------------------------
! CALCULATE EVAPORATION
!-----------------------------------------------------------------------
!
! DEPENDS ON: evp
        CALL EVP (NPNTS,RAIN,T,CCA,RHO,DELQ,DELP,EVAP_RAIN,             &
     &            BEVAP,1,CLDAREA,PK)
!
! DEPENDS ON: evp
        CALL EVP (NPNTS,SNOW,T,CCA,RHO,DELQ,DELP,SUB_SNOW,              &
     &            BEVAP,2,CLDAREA,PK)
!
!-----------------------------------------------------------------------
! CALCULATE TEMPERATURE AND MIXING RATIO IF LAYER BROUGHT TO
! SATURATION BY EVAPORATION AND SUBLIMATION
!-----------------------------------------------------------------------
!
! DEPENDS ON: satcal
      CALL SATCAL(NPNTS,T,TH,PK,QS,THS,K,EXK,Q,TH)
!
!
      DO I=1,NPNTS
        IF (BEVAP(I)) THEN
          DTHBYDT_EVP(I) = -((LC*EVAP_RAIN(I))+((LC+LF)*SUB_SNOW(I)))*G/&
     &                   (CP*EXK(I)*DELP(I))
          DQBYDT_EVP(I) = (EVAP_RAIN(I)+SUB_SNOW(I))*G/DELP(I)
!
          DTHBYDT_SAT(I) = (THS(I)-TH(I))/TIMESTEP
!
          IF (DTHBYDT_EVP(I) <  DTHBYDT_SAT(I)) THEN
!
!---------------------------------------------------------------------
!  ADJUST EVAPORATION AND SUBLIMATION RATES TO GIVE SATURATION
!---------------------------------------------------------------------
!
            FACTOR(I) = DTHBYDT_SAT(I)/DTHBYDT_EVP(I)
            DTHBYDT_EVP(I) = DTHBYDT_SAT(I)
            DQBYDT_EVP(I) = DQBYDT_EVP(I)*FACTOR(I)
            EVAP_RAIN(I) = EVAP_RAIN(I)*FACTOR(I)
            SUB_SNOW(I) = SUB_SNOW(I)*FACTOR(I)
          END IF
!
!---------------------------------------------------------------------
!  UPDATE INCREMENTS AND RAINFALL AND ADJUST BACK TO GRIDBOX MEANS
!---------------------------------------------------------------------
!
          DTHBYDT(I) = DTHBYDT(I)+DTHBYDT_EVP(I)*CCA(I)*CLDAREA
          DQBYDT(I) = DQBYDT(I)+DQBYDT_EVP(I)*CCA(I)*CLDAREA
          RAIN(I) = RAIN(I)-EVAP_RAIN(I)*CCA(I)*CLDAREA
          SNOW(I) = SNOW(I)-SUB_SNOW(I)*CCA(I)*CLDAREA
        END IF
      END DO
!
      RETURN
      END SUBROUTINE PEVP_BCB
!
#endif
