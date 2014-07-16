
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE EVP----------------------------------------------------
!LL
!LL  PURPOSE : CALCULATES THE EVAPORATION OF PRECIPITATION
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL  4.0    5/05/95   New deck at version 4.0 to include pressure
!LL                   dependency into calculation of evaporation of
!LL                   convective precipitation.
!LL                   Pete Inness.
!LL  4.2    Oct. 96   T3E migration: *DEF CRAY removed, HF functions
!LL                       replaced.
!LL                                  S.J.Swarbrick
!LL  4.3    Feb. 97   T3E optimisation of powers & sqrt
!LL                                  D.Salmond & S.J.Swarbrick
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!LL  5.3   7/8/01     Prevent negative evaporation of falling snow.  WJI
!    5.5  17/04/03  Removal of reference to obsolete section
!                   A05_3B. T.White
!    6.2  03/02/05  Added section 5A. R A Stratton
!    6.4  04/12/06  Removed old T3E optimisation to make things
!                   more vector friendly. P.Selwood
!    6.4  12/12/06  Removed def 3C. R A Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL  VERSION NO. 4  DATED 23/7/92
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
      SUBROUTINE EVP(NPNTS,PRECIP,TEVP,CCA,RHO,DELQ,DELPKM1,EVAP,       &
     &               BEVAP,IPHASE,AREA_FAC,PKM1)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS USED IN THIS SUBROUTINE
!-----------------------------------------------------------------------
!
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! DDEVAP start
      ! exponents used in calculation of evaporation of liquid
      REAL,PARAMETER :: P_LQ1 = 0.52
      REAL,PARAMETER :: P_LQ2 = 0.67

      ! exponents used in calculation of evaporation of ice
      REAL,PARAMETER :: P_ICE1 = 0.55
      REAL,PARAMETER :: P_ICE2 = 0.76

      ! exponents and constants associated with density term in
      ! evaporation of liquid
      REAL,PARAMETER :: RHO_LQP1 = 0.26
      REAL,PARAMETER :: RHO_LQP2 = 0.59
      REAL,PARAMETER :: RHO_LQA  = 108.80
      REAL,PARAMETER :: RHO_LQB  = 830.73

      ! exponents and constants associated with density term in
      ! evaporation of ice
      REAL,PARAMETER :: RHO_ICP1 = 0.28
      REAL,PARAMETER :: RHO_ICP2 = 0.63
      REAL,PARAMETER :: RHO_ICEA = 1569.52
      REAL,PARAMETER :: RHO_ICEB = 32069.02
! DDEVAP end
! DDEVPLQ start

      ! constants used in quadratic formula for evaporation of liquid
      REAL,PARAMETER:: LQ_A = 2.008E-9
      REAL,PARAMETER:: LQ_B = -1.385E-6
      REAL,PARAMETER:: LQ_C = 2.424E-4

! DDEVPLQ end
! DDEVPICE start

      ! constants used in quadratic formula for evaporation of ice
      REAL,PARAMETER:: ICE_A = -5.2E-9
      REAL,PARAMETER:: ICE_B = 2.5332E-6
      REAL,PARAMETER:: ICE_C = -2.911E-4

! DDEVPICE end
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
      INTEGER I                 ! LOOP COUNTER
!
      INTEGER NPNTS             ! IN VECTOR LENGTH
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
!
      REAL DELQ(NPNTS)          ! IN CHANGE IN HUMIDITY MIXING
                                !    RATIO ACROSS LAYER K (KG/KG)
!
      REAL TEVP(NPNTS)          ! IN TEMPERATURE OF LAYER K (K)
!
      LOGICAL BEVAP(NPNTS)      ! IN MASK FOR POINTS WHERE EVAPORATION
                                !    TAKES PLACE
!
      REAL PRECIP(NPNTS)        ! IN AMOUNT OF PRECIPITATION(KG/M**2/S)
!
      REAL DELPKM1(NPNTS)       ! IN CHANGE IN PRESSURE ACROSS
                                !    LAYER K-1
!
      REAL CCA(NPNTS)           ! IN CONVECTIVE CLOUD AMOUNT
!
      REAL RHO(NPNTS)           ! IN DENSITY OF AIR
!
      INTEGER IPHASE            ! IN INDICATION FOR RAIN (1), OR
                                !    SNOW (2)
!
      REAL PKM1(NPNTS)          ! IN PRESSURE AT LEVEL KM1
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!-----------------------------------------------------------------------
!
      REAL EVAP(NPNTS)   ! OUT EVAPORATION
!
!-----------------------------------------------------------------------
! EXTERNAL ROUTINES
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE LOCALLY DEFINED
!-----------------------------------------------------------------------
!
      REAL ECON          ! QUADRATIC TERM
!
      REAL C1            ! CONSTANT
!
      REAL C2            ! CONSTANT
!
      REAL LRATE         ! LOCAL RATE OF PRECIPITATION
!
      REAL CA            ! LOCAL CLOUD AREA
!
      REAL AREA_FAC      ! FRACTION OF CONVECTIVE CLOUD AMOUNT TO GIVE
                         ! LOCAL CLOUD AREA
      real tl1,ti1
!
!-----------------------------------------------------------------------
! START OF ROUTINE
!-----------------------------------------------------------------------
!
      tl1=0.5*P_LQ1
      ti1=0.5*P_ICE1
!
      IF (IPHASE == 1) THEN        ! RAIN
!
      DO I=1,NPNTS
       IF (BEVAP(I)) THEN
         IF (PRECIP(I)  >   0.0) THEN
           ECON = ((LQ_A*TEVP(I)+LQ_B)*TEVP(I)+LQ_C)*                   &
     &                         (100000.0/PKM1(I))
           CA = AREA_FAC*CCA(I)
           LRATE = PRECIP(I)/CA
           C1 = RHO_LQA*CA*(LRATE*LRATE*RHO(I))**tl1
           C2 = RHO_LQB*CA*LRATE**P_LQ2*RHO(I)**RHO_LQP2
           EVAP(I) = MIN(ECON*(C1+C2)*DELQ(I)*DELPKM1(I)/G,LRATE)
         ELSE
           EVAP(I) = 0.0
         END IF
       END IF
      END DO
!
      ELSE IF (IPHASE == 2) THEN        ! SNOW
!
!
      DO I=1,NPNTS
       IF (BEVAP(I)) THEN
         IF (PRECIP(I)  >   0.0) THEN
           IF(TEVP(I) <= 243.58) THEN
             ECON = 1.7405E-5*(100000.0/PKM1(I))
           ELSE
             ECON = ((ICE_A*TEVP(I)+ICE_B)*TEVP(I)+ICE_C)*              &
     &                    (100000.0/PKM1(I))
           END IF
           CA = AREA_FAC*CCA(I)
           LRATE = PRECIP(I)/CA
           C1 = RHO_ICEA*CA*(LRATE*LRATE*RHO(I))**ti1
           C2 = RHO_ICEB*CA*LRATE**P_ICE2*RHO(I)**RHO_ICP2
           EVAP(I)=MAX(0.,MIN(ECON*(C1+C2)*DELQ(I)*DELPKM1(I)/G,LRATE))

         ELSE
           EVAP(I) = 0.0
         END IF
       END IF
      END DO
!
      ENDIF
!
      RETURN
      END SUBROUTINE EVP
!
