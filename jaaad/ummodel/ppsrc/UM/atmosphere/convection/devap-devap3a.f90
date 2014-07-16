
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DEVAP--------------------------------------------------
!LL
!LL  PURPOSE : EVAPORATION ROUTINE
!LL
!LL            CARRIES OUT EVAPORATION AND UPDATES PRECIPITATION
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL  4.0    5/05/95   New deck at version 4.0 to include pressure
!LL                   dependency into calculation of evaporation of
!LL                   convective precipitation, and to introduce
!LL                   traps for negative precipitation.
!LL                   Pete Inness.
!LL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!    5.5   17/04/03  Removal of reference to obsolete section
!                    A05_3B. T.White
!    6.2  03/02/05   Added section 5A. R A Stratton
!    6.4  12/12/06   Removed def 3C. R A Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL  VERSION NO. 4  DATED 5/2/92
!LL
!LL  SYSTEM TASK : P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE DEVAP(NPNTS,THDD_K,THDD_KM1,QDD_KM1,THDDS,QDDS,        &
     &                 FLX_DD_KM1,EXK,EXKM1,QSATDD,RAIN,SNOW,           &
     &                 DELPKM1,BDDWT_KM1,CCA,PKM1)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS
!-----------------------------------------------------------------------
!
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
! DDAREA start

      REAL,PARAMETER:: DDCLDFRA=0.5 ! fractional cloud area of dd

! DDAREA end
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
!
      INTEGER I               ! LOOP COUNTER
!
      INTEGER NPNTS           ! IN VECTOR LENGTH
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
      REAL THDDS(NPNTS)       ! IN SATURATED POTENTIAL
                              !    TEMPERATURE OF DOWNDRAUGHT
                              !    (K)
!
      REAL QDDS(NPNTS)        ! IN MIXING RATIO OF SATURATED
                              !    DOWNDRAUGHT (KG/KG)
!
      REAL FLX_DD_KM1(NPNTS)  ! IN DOWNDRAUGHT MASS FLUX IN
                              !    LAYER K-1 (PA/S)
!
      REAL THDD_K(NPNTS)      ! IN POTENTIAL TEMPERATURE OF
                              !    DOWNDRAUGHT IN LAYER K (K)
!
      REAL EXK(NPNTS)         ! IN EXNER RATIO OF LAYER K
!
      REAL EXKM1(NPNTS)       ! IN EXNER RATIO OF LAYER K-1
!
      REAL QSATDD(NPNTS)      ! IN SATURATED DOWNDRAUGHT
                              !    MIXING RATIO (KG/KG)
!
      REAL DELPKM1(NPNTS)     ! IN CHANGE IN PRESSURE ACROSS
                              !    LAYER K-1 (PA)
!
      LOGICAL BDDWT_KM1(NPNTS)! IN MASK WHERE PRECIPITATION IN
                              !    DOWNDRAUGHT IS LIQUID
!
      REAL CCA(NPNTS)         ! IN CONVECTIVE CLOUD AMOUNT
!
      REAL PKM1(NPNTS)        ! IN PRESSURE OF LAYER K-1
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!-----------------------------------------------------------------------
!
      REAL THDD_KM1(NPNTS)    ! INOUT
                              ! IN  POTENTIAL TEMPERATURE OF
                              !     DOWNDRAUGHT IN LAYER K-1 (K)
                              ! OUT UPDATED POTENTIAL TEMPERATURE
                              !     OF DOWNDRAUGHT IN LAYER K-1
                              !     AFTER EVAPORATION OR
                              !     SATURATION (K)
!
      REAL QDD_KM1(NPNTS)     ! INOUT
                              ! IN  MODEL MIXING RATIO OF
                              !     DOWNDRAUGHT IN LAYER K-1
                              !     (KG/KG)
                              ! OUT UPDATED MODEL MIXING RATIO
                              !     OF DOWNDRAUGHT IN LAYER K-1
                              !     AFTER EVAPORATION OR
                              !     SATURATION (KG/KG)
!
      REAL RAIN(NPNTS)        ! INOUT
                              ! IN  AMOUNT OF RAIN (KG/M**2/S)
                              ! OUT UPDATED RAINFALL (KG/M**2/S)
!
      REAL SNOW(NPNTS)        ! INOUT
                              ! IN  AMOUNT OF SNOW (KG/M**2/S)
                              ! OUT UPDATED SNOWFALL (KG/M**2/S)
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE LOCALLY DEFINED
!-----------------------------------------------------------------------
!
!
      REAL TEVP(NPNTS)        ! TEMPERATURE USED IN EVAPORATION
                              ! CALCULATION (K)
!
      LOGICAL BEVAP(NPNTS)    ! MASK FOR THOSE POINTS AT WHICH
                              ! EVAPORATION CALCULATION IS TO
                              ! BE CARRIED OUT
!
      LOGICAL BSAT(NPNTS)     ! MASK FOR THOSE POINTS WHICH
                              ! ARE SUBSATURATED
!
      REAL EVAP_RAIN(NPNTS)   ! AMOUNT OF EVAPORATION OF RAIN
!
      REAL SUB_SNOW(NPNTS)    ! AMOUNT OF SNOW SUBLIMATION
!
      REAL DELQ(NPNTS)        ! DIFFERENCE IN MIXING RATIOS
                              ! (KG/KG)
!
      REAL DELTH(NPNTS)       ! INCREMENT TO DOWNDRAUGHT POTENTIAL
                              ! TEMPERATURE IN LAYER K-1 DUE TO
                              ! EVAPORATION
!
      REAL DELQE(NPNTS)       ! INCREMENT TO DOWNDRAUGHT MIXING RATIO
                              ! IN LAYER K-1 DUE TO EVAPORATION
!
      REAL DELTHS(NPNTS)      ! SATURATED POTENTIAL TEMPERATURE MINUS
                              ! POTENTIAL TEMPERATURE OF DOWNDRAUGHT
!
      REAL FACTOR(NPNTS)      ! DELTHS / DELTH
!
      REAL PINCR(NPNTS)       ! INCREASE IN PRECIPITATION IF PARCEL
                              ! SUPERSATURATES
!
      REAL RHO(NPNTS)         ! DENSITY OF AIR IN PARCEL
!
!
!-----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!-----------------------------------------------------------------------
!
      EXTERNAL EVP
!
!-----------------------------------------------------------------------
! CHECK IF EVAPORATION POSSIBLE
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
       DELQ(I) = QSATDD(I)-QDD_KM1(I)
!
       BEVAP(I) =((RAIN(I) >  0.0) .OR. (SNOW(I) >  0.0))               &
     &             .AND. (DELQ(I) >  0.0)
       BSAT(I) = DELQ(I)  <   0.0
!
!-----------------------------------------------------------------------
! CALCULATE TEMPERATURE USED IN CALCULATION OF EVAPORATION CONSTANTS
! BASED ON TEMPERATURE OF PARCEL AFTER UNSATURATED DESCENT
!-----------------------------------------------------------------------
!
        IF (BEVAP(I)) THEN
          TEVP(I) = ((THDD_K(I)*EXK(I))+(THDD_KM1(I)*EXKM1(I)))*0.5
          RHO(I) = PKM1(I) / (R*TEVP(I))
        END IF
      END DO
!
!-----------------------------------------------------------------------
! EVAPORATION CALCULATION - CALCULATE RATES FOR RAIN AND SNOW
!-----------------------------------------------------------------------
!
! DEPENDS ON: evp
      CALL EVP(NPNTS,RAIN,TEVP,CCA,RHO,DELQ,DELPKM1,EVAP_RAIN,          &
     &         BEVAP,1,DDCLDFRA,PKM1)
!
! DEPENDS ON: evp
      CALL EVP(NPNTS,SNOW,TEVP,CCA,RHO,DELQ,DELPKM1,SUB_SNOW,           &
     &         BEVAP,2,DDCLDFRA,PKM1)
!
      DO I=1,NPNTS
       IF (BEVAP(I)) THEN
!
!-----------------------------------------------------------------------
! ADJUST EVAPORATION AND SUBLIMATION RATES BACK TO GRID BOX MEANS
!-----------------------------------------------------------------------
!
       EVAP_RAIN(I) = EVAP_RAIN(I) * CCA(I) * DDCLDFRA
       SUB_SNOW(I) = SUB_SNOW(I) * CCA(I) * DDCLDFRA
!
!-----------------------------------------------------------------------
! CHECK IF PARCEL SUPERSATURATED
!-----------------------------------------------------------------------
!
        DELTH(I) = -((LC*EVAP_RAIN(I))+((LC+LF)*SUB_SNOW(I)))*G/        &
     &           (CP*EXKM1(I)*FLX_DD_KM1(I))
        DELQE(I) = (EVAP_RAIN(I)+SUB_SNOW(I))*G/FLX_DD_KM1(I)
!
        DELTHS(I) = THDDS(I)-THDD_KM1(I)
        IF (DELTH(I) <  DELTHS(I)) THEN
!
!-----------------------------------------------------------------------
! ADJUST EVAP AND SUBLIMATION RATES TO GIVE SATURATION
!-----------------------------------------------------------------------
!
          FACTOR(I) = DELTHS(I)/DELTH(I)
          DELTH(I) = DELTHS(I)
          DELQE(I) = DELQE(I)*FACTOR(I)
          EVAP_RAIN(I) = EVAP_RAIN(I)*FACTOR(I)
          SUB_SNOW(I) = SUB_SNOW(I)*FACTOR(I)
        END IF
!
!-----------------------------------------------------------------------
! UPDATE T,Q AND PRECIPITATION
!-----------------------------------------------------------------------
!
        RAIN(I) = RAIN(I)-EVAP_RAIN(I)
        IF (RAIN(I) <  0.0) RAIN(I)=0.0
        SNOW(I) = SNOW(I)-SUB_SNOW(I)
        IF (SNOW(I) <  0.0) SNOW(I)=0.0
        THDD_KM1(I) = THDD_KM1(I)+DELTH(I)
        QDD_KM1(I) = QDD_KM1(I)+DELQE(I)
!
!-----------------------------------------------------------------------
! PARCEL IS SUPERSATURATED BEFORE EVAPORATION OCCURS
! BRING PARCEL TO SATURATION AND PRECIPITATE WATER
!-----------------------------------------------------------------------
!
      ELSE IF (BSAT(I)) THEN
         PINCR(I) = (QDD_KM1(I)-QDDS(I))*FLX_DD_KM1(I)/G
         QDD_KM1(I) = QDDS(I)
         THDD_KM1(I) = THDDS(I)
         IF (BDDWT_KM1(I)) THEN
           RAIN(I) = RAIN(I)+PINCR(I)
         ELSE
           SNOW(I) = SNOW(I)+PINCR(I)
         END IF
      END IF
      END DO
!
      RETURN
      END SUBROUTINE DEVAP
!
