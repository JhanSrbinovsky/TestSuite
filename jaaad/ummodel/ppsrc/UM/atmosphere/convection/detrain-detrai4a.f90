
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DETRAIN------------------------------------------------
!LL
!LL  PURPOSE : FORCED DETRAINMENT CALCULATION
!LL
!LL            SUBROUTINE THP_DET CALCULATES THE POTENTIAL
!LL            TEMPERATURE OF THE PARCEL IN LAYER K+1
!LL            AFTER FORCED DETRAINMENT
!LL
!LL            SUBROUTINE THETAR CALCULATES THE POTENTIAL TEMPERATURE
!LL            OF THE AIR IN LAYER K UNDERGOING FORCED DETRAINMENT
!LL
!LL            SUBROUTINE DET_RATE CALCULATES THE FORCED DETRAINMENT
!LL            RATE OF THE ENSEMBLE IN LAYER K
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL
!LL   5.4    6/8/02    New deck created for version 4A of convection
!LL                    scheme.            Gill Martin
!     6.2   03/02/05  Added section 5A R.A.Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
!LL  VERSION NO. 1
!LL
!LL  LOGICAL COMPONENTS COVERED: P27
!LL
!LL  SYSTEM TASK :
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE DETRAIN (NPNTS,THEK,QEK,THPK,QPK,QSEK,DQSK,BGMK,       &
     &                     THEKP1,QEKP1,THPKP1,QPKP1,QSEKP1,DQSKP1,     &
     &                     BGMKP1,BWKP1,XSQKP1,                         &
     &                     DELTAK,THRK,QRK,EKP14,EKP34,PK,PKP1,         &
     &                     XSBMIN,                                      &
     &                     EXK,EXKP1)
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NPNTS          ! IN VECTOR LENGTH
!
      INTEGER I              ! LOOP COUNTER
!
      INTEGER NREDO          ! NUMBER OF POINTS FOR WHICH FORCED
                             ! DETRAINMENT CALCULATION MUST BE
                             ! AS THE PROCESSES EITHER CAUSES THE
                             ! PARCEL TO BECOME SATURATED OR
                             ! SUB-SATURATED
!
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      REAL THEK(NPNTS)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K (K)
!
      REAL THEKP1(NPNTS)     ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
!
      REAL QEK(NPNTS)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL QEKP1(NPNTS)      ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL DQSKP1(NPNTS)     ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT IN LAYER K+1
                             !    (KG/KG/K)
!
      REAL THPK(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K (K)
!
      REAL QPK(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
!
      REAL QSEK(NPNTS)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL DQSK(NPNTS)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT OF LAYER K
                             !    (KG/KG/K)
!
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
!
      LOGICAL BGMK(NPNTS)    ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K
!
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL EKP34(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+3/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL EXKP1(NPNTS)      ! IN EXNER RATIO AT LEVEL K+1
!
      REAL EXK(NPNTS)        ! IN EXNER RATIO AT LEVEL K
!
      REAL PKP1(NPNTS)       ! IN PRESSURE AT LEVEL K+1 (PA)
!
      REAL PK(NPNTS)         ! IN PRESSURE AT LEVEL K (PA)
!
      REAL XSBMIN(NPNTS)     ! IN THRESHOLD FOR FORCED DETRAINMENT
!                            !    Function of delta P
!
!
!-----------------------------------------------------------------------
! VARIABLES WHICH INPUT AND OUTPUT
!-----------------------------------------------------------------------
!
      REAL THPKP1(NPNTS)     ! INOUT
                             ! IN  PARCEL POTENTIAL TEMPERATURE IN
                             !     LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (K)
                             ! OUT ADJUSTED PARCEL POTENTIAL
                             !     IN LAYER K+1 AFTER FORCED
                             !     DETRAINMENT (K)
!
      REAL QPKP1(NPNTS)      ! INOUT
                             ! IN  PARCEL MIXING RATIO IN
                             !     LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (KG/KG)
                             ! OUT ADJUSTED PARCEL POTENTIAL
                             !     IN LAYER K+1 AFTER FORCED
                             !     DETRAINMENT (KG/KG)
!
      REAL XSQKP1(NPNTS)     ! INOUT
                             ! IN  EXCESS WATER IN PARCEL AFTER
                             !     LIFTING FROM LAYER K TO K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING
                             !     (KG/KG)
                             ! OUT EXCESS WATER IN PARCEL IN LAYER
                             !     K+1 AFTER FORCED DETRAINMENT
                             !     (KG/KG)
!
      LOGICAL BGMKP1(NPNTS)  ! INOUT
                             ! IN  MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING
                             ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 AFTER
                             !     FORCED DETRAINMENT
!
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!-----------------------------------------------------------------------
!
      REAL THRK(NPNTS)       ! OUT PARCEL DETRAINMENT POTENTIAL
                             !     TEMPERATURE IN LAYER K (K)
!
      REAL QRK(NPNTS)        ! OUT PARCEL DETRAINMENT MIXING RATIO
                             !     IN LAYER K (KG/KG)
!
      REAL DELTAK(NPNTS)     ! OUT PARCEL FORCED DETRAINMENT RATE
                             !     IN LAYER K
!
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!
      LOGICAL BDETK(NPNTS)   ! MASK FOR PARCELS WHICH ARE
                             ! UNDERGOING FORCED DETRAINMENT
                             ! IN THEIR ASCENT FROM LAYER K
                             ! TO K+1
!
      REAL XSQR(NPNTS)       ! EXCESS PARCEL WATER VAPOUR
                             ! DURING DETRAINMENT (KG/KG)
!
      REAL THPKP1W(NPNTS) ,                                             &
                             ! TEMPORARY STOREAGE FOR PARCEL
     &     QPKP1W(NPNTS) ,                                              &
                             ! POTENTIAL TEMPERATURE (K), MIXING
     &     XSQK1W(NPNTS)     ! RATIO (KG/KG), EXCESS WATER VAPOUR
      LOGICAL BGKP1W(NPNTS)  ! (KG/KG) AND MASK FOR SATURATION
                             ! IN LAYER K+1
!
      LOGICAL BRECAL(NPNTS)  ! MASK FOR THOSE POINTS AT WHICH THE
                             ! THE DETRAINMENT CALCULATION NEEDS
                             ! REPEATING
!
      REAL TT(NPNTS)         ! TEMPORARY STORE FOR TEMPERATURE
                             ! FOR THE CALCULATION OF SATURATED
                             ! MIXING RATIO (K)
!
      REAL EPSS              ! (1+EKP14)*(1+EKP34)
!
!----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!----------------------------------------------------------------------
!
      EXTERNAL THP_DET,QSAT,THETAR,DET_RATE
!
!*---------------------------------------------------------------------
!
      DO 10 I=1,NPNTS
!
!----------------------------------------------------------------------
! AT START OF ROUTINE FORCED DETARINMENT DONE AT ALL POINTS SO
! SET ARRAY BDETK EQUAL TO .TRUE.
! SET FORCED DETRAINMENT RATE EQUAL TO ZERO
!----------------------------------------------------------------------
!
       BDETK(I) = .TRUE.
       DELTAK(I) = 0.0
!
!-----------------------------------------------------------------------
!   SAVE THE CURRENT VALUES OF QPKP1, XSQKP1 AND BGMKP1
!-----------------------------------------------------------------------
!
       THPKP1W(I) = THPKP1(I)
       QPKP1W(I) = QPKP1(I)
       XSQK1W(I) = XSQKP1(I)
       BGKP1W(I) = BGMKP1(I)
!
!-----------------------------------------------------------------------
!   ADD THE EXCESS WATER VAPOUR BACK INTO THE DETRAINING PARCELS
!-----------------------------------------------------------------------
!
       QPKP1(I) = QPKP1(I) + XSQKP1(I)
   10 CONTINUE
!L
!L----------------------------------------------------------------------
!L  CALCULATE THE ENSEMBLE AVERAGE POTENTIAL TEMPERATURE IN LAYER K+1
!L  AT THE POINTS WHERE DETRAINMENT IS TAKING PLACE
!L
!L  SUBROUTINE THP_DET
!L
!L  UM DOCUMENTATION PAPER P27
!L  SECTION (6), EQUATION (28)
!L----------------------------------------------------------------------
!L
! DEPENDS ON: thp_det
      CALL THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,DQSKP1,      &
     &              BGMKP1,BDETK,PK,PKP1,XSBMIN)
!L
!L---------------------------------------------------------------------
!L  CHECK TO SEE IF SUFFICIENT EXCESS WATER VAPOUR IN THE
!L  INITIAL DRY ASCENT TO ALLOW PARCEL TO BE SATURATED
!L  IN LAYER K+1 AFTER FORCED DETRAINMENT
!L
!L  UM DOCUMENTATION PAPER P27
!L  SECTION (6), EQUATION (29)
!L
!L  NOTE : ONLY ALLOW PARCEL TO BE SATURATED IN LAYER K+1 IF
!L         SATURATED INITIALLY.  IT IS POSSIBLE FOR SMALL
!L         SUPERSATURATIONS TO IF SUBROUTINE LATENT_H CAUSES
!L         PARCEL TO BE COME UNSATURATED.  IN THIS CASE TREAT
!L         THE PARCEL AS UNSATURATED IN LAYER K+1
!L---------------------------------------------------------------------
!L
!
!-----------------------------------------------------------------------
!   CALCULATE THE EXCESS WATER VAPOUR IN LAYER K+1 AND RECALCULATE
!   BGMKP1 AND QPKP1.
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
! CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
! PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
! MIXING RATIO
!-----------------------------------------------------------------------
!
      DO 25 I = 1,NPNTS
       TT(I) = THPKP1(I)*EXKP1(I)
   25 CONTINUE
! DEPENDS ON: qsat
      CALL QSAT (XSQKP1,TT,PKP1,NPNTS)
!
      DO 30 I=1,NPNTS
       XSQKP1(I) = QPKP1(I) - XSQKP1(I)
!
       BRECAL(I) = BGMKP1(I)
!
!----------------------------------------------------------------------
! ONLY ALLOW PARCEL TO BE SATURATED IN INITIAL BGMKP1 = .TRUE.
! (STORED IN BRECAL AT THIS POINT)
!----------------------------------------------------------------------
!
       IF ( BGMK(I) .OR.( (XSQKP1(I)  >   0.) .AND. BRECAL(I) ) ) THEN
         BGMKP1(I) = .TRUE.
       ELSE
         BGMKP1(I) = .FALSE.
         XSQKP1(I) = 0.0
       END IF
!
       QPKP1(I) = QPKP1(I) - XSQKP1(I)
!L
!L----------------------------------------------------------------------
!L  RECALCULATE THE ENSEMBLE AVERAGE POTENTIAL TEMPERATURE AT POINTS
!L  WHERE THE ENSEMBLE HAS BECOME UNSATURATED.
!L
!L  UM DOCUMENTATION PAPER P27
!L  SECTION (6), EQUATION (28)
!L----------------------------------------------------------------------
!L
       BRECAL(I) = BDETK(I) .AND. BRECAL(I) .AND. .NOT.BGMKP1(I)
   30 CONTINUE
!
! DEPENDS ON: thp_det
      CALL THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,DQSKP1,      &
     &             BGMKP1,BRECAL,PK,PKP1,XSBMIN)
!L
!L----------------------------------------------------------------------
!L  BECAUSE OF THE REMOVAL OF LATENT HEATING, THE NEW PARCEL POTENTIAL
!L  TEMPERATURE MAY BE LOWER THAN ITS VALUE BEFORE THE DETRAINMENT
!L  CALCULATION. IN THIS CASE ABANDON THE DETRAINMENT CALCULATION.
!L----------------------------------------------------------------------
!L
      DO 90 I=1,NPNTS
       BDETK(I) = THPKP1(I)  >   THPKP1W(I)
   90 CONTINUE
!L
!L----------------------------------------------------------------------
!L  CALCULATE THE POTENTIAL TEMPERATURE AND MIXING RATIO  OF DETRAINING
!L  AIR AND THE EXCESS WATER VAPOUR CONDESED FROM DETRAINING AIR
!L
!L  UM DOCUMENTATION PAPER P27
!L  SECTION (6), EQUATION (26)
!L----------------------------------------------------------------------
!L
! DEPENDS ON: thetar
      CALL THETAR (NPNTS,THRK,QRK,XSQR,BGMK,THEK,QEK,QPK,QSEK,DQSK,     &
     &             BWKP1,EXK,PK)
!L
!L----------------------------------------------------------------------
!L  CALCULATE THE DETRAINMENT RATE, DELTAK.
!L
!L  UM DOCUMENTATION PAPER P27
!L  SECTION (6), EQUATION (31)
!L----------------------------------------------------------------------
!L
! DEPENDS ON: det_rate
      CALL DET_RATE (NPNTS,DELTAK,THRK,XSQR,THPK,THEK,THEKP1,           &
     &             XSQKP1,THPKP1,BWKP1,BDETK,EKP14,EKP34,EXK,EXKP1)
!
      NREDO = 0
!L
!L----------------------------------------------------------------------
!L  ADD WATER VAPOUR WHICH WAS REMOVED FROM DETRAINING AIR INTO XSQKP1
!L
!L  UM DOCUMENTATION PAPER P27
!L  SECTION 86), EQUATION (11C)
!L----------------------------------------------------------------------
!L
      DO 120 I=1,NPNTS
!
       EPSS = (1.+EKP14(I))*(1.+EKP34(I))
!
       IF (BDETK(I))                                                    &
     & XSQKP1(I) = XSQKP1(I) + (DELTAK(I)*XSQR(I)/                      &
     &               (EPSS*(1.-DELTAK(I))))
!L
!L----------------------------------------------------------------------
!L  IF THE EXCESS WATER VAPOUR IN LAYER K+1 IS LESS THAN ZERO
!L  I.E. THE PARCEL HAS BECOME UNSATURATED THROUGH THE FORCED
!L  DETRAINMENT PROCESS THEN ABANDON THE CALCULATION
!L----------------------------------------------------------------------
!L
       BRECAL(I) = BGMKP1(I)
!
       BGMKP1(I) = XSQKP1(I)  >   0.
!
       BRECAL(I) = BDETK(I) .AND. BRECAL(I) .AND. .NOT.BGMKP1(I)
!
       IF (BRECAL(I)) THEN
          QPKP1(I)  = QPKP1(I) + XSQKP1(I)                              &
     &               - (DELTAK(I)*XSQR(I)/(EPSS*(1.-DELTAK(I))))
          XSQKP1(I) = 0.
       ENDIF
!
!----------------------------------------------------------------------
! COUNT POINTS AT WHICH DETRAINMENT CALCULATION NEEDS REPEATING
!----------------------------------------------------------------------
!
       IF (BRECAL(I)) NREDO = NREDO + 1
  120 CONTINUE
!L
!L---------------------------------------------------------------------
!L  REPEAT CALCULATION OF PARCEL POTENTIAL TEMPERATURE, DETRAINMENT
!L  RATE AND EXCESS PARCEL WATER IF THE PARCEL BECOMES UNSATURATED
!L  IN LAYER K+1 AFTER FORCED DETARINMENT
!L---------------------------------------------------------------------
!L
      IF (NREDO  >   0) THEN
!
!----------------------------------------------------------------------
!  CALCULATE NEW PARCEL POTENTIAL TEMPERATURE IN LAYER K+1
!  AFTER FORCED DETRAINMENT
!----------------------------------------------------------------------
!
! DEPENDS ON: thp_det
        CALL THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,DQSKP1,    &
     &               BGMKP1,BRECAL,PK,PKP1,XSBMIN)
!
!----------------------------------------------------------------------
!  CHECK IF FORCED DETRAINMENT STILL POSSIBLE AND RESET RECALCUATION
!  MASK TO FALSE IF IT IS NOT
!----------------------------------------------------------------------
!
        DO 130 I=1,NPNTS
          IF (BRECAL(I)) THEN
            BDETK(I) = THPKP1(I)  >   THPKP1W(I)
            BRECAL(I) = BDETK(I)
          END IF
  130   CONTINUE
!
!----------------------------------------------------------------------
!  RCALCULATE FORCED DETRAINEMNT RATE
!----------------------------------------------------------------------
!
! DEPENDS ON: det_rate
        CALL DET_RATE (NPNTS,DELTAK,THRK,XSQR,THPK,THEK,THEKP1,         &
     &             XSQKP1,THPKP1,BWKP1,BRECAL,EKP14,EKP34,EXK,EXKP1)
!
!----------------------------------------------------------------------
!  RECALCULATE EXCESS WATER VAPOUR IN LAYER K+1
!  AFTER FORCED DETRAINMENT
!----------------------------------------------------------------------
!
        DO 140 I=1,NPNTS
         IF (BRECAL(I)) THEN
            EPSS = (1.+EKP14(I))*(1.+EKP34(I))
            XSQKP1(I) = XSQKP1(I) + (DELTAK(I)*XSQR(I)/                 &
     &                       (EPSS*(1.-DELTAK(I))))
         END IF
  140   CONTINUE
!
      END IF
!L
!L----------------------------------------------------------------------
!L  MAKE SURE THAT THE DETRAINMENT RATE IS BETWEEN 0 AND 1
!L
!L  IF <0 THEN NO DETRAINMENT OCCURS AND ORIGINAL VALUES ARE
!L  RESTORED
!L
!L  IF >1 THEN SET TO 1 AND THRK = THPK, QRK = QPK AND VALUES
!L  IN LAYER K+1 ARE RESTORED.  ALTHOUGH THESE ARE NOT USED
!L  IN ANY THERMODYNAMIC CALCULATION THEY ARE USED TO SPECIFY
!L CLOUD TOP IN SUBROUTIBE CONRAD
!L----------------------------------------------------------------------
!L
      DO 180 I=1,NPNTS
!
       IF (BDETK(I)) THEN
!
        IF (DELTAK(I) <= 0.0) THEN
           BDETK(I) = .FALSE.
           THPKP1 (I) = THPKP1W(I)
           QPKP1 (I) = QPKP1W(I)
           XSQKP1(I) = XSQK1W(I)
           BGMKP1(I) = BGKP1W(I)
           DELTAK(I) = 0.0
        ELSE IF (DELTAK(I) >  1.0) THEN
           DELTAK(I) = 1.0
           THRK(I) = THPK(I)
           QRK(I) = QPK(I)
           THPKP1 (I) = THPKP1W(I)
           QPKP1 (I) = QPKP1W(I)
           XSQKP1(I) = XSQK1W(I)
           BGMKP1(I) = BGKP1W(I)
        END IF
!
       ELSE

        THPKP1 (I) = THPKP1W(I)
        QPKP1 (I) = QPKP1W(I)
        XSQKP1(I) = XSQK1W(I)
        BGMKP1(I) = BGKP1W(I)
        DELTAK(I) = 0.0

       ENDIF
  180  CONTINUE
!
      RETURN
      END SUBROUTINE DETRAIN
