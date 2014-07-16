
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE THETAR-------------------------------------------------
!LL
!LL  PURPOSE : CALCULATES THE POTENTIAL TEMPERATURE OF THE DETRAINING
!LL            AIR IN LAYER K AND ALSO THE DIFFERENCE IN THE
!LL            WATER VAPOUR CONTENT OF THE DETRAINING AIR FROM THAT
!LL            OF THE MEAN PARCEL IN LAYER K
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!     5.5  17/04/03   Removal of references to obsolete sections
!                     A05_2A,2C,3B. T.White
!LL
!LL  4.1  6/6/96     Extra check added to ensure that
!LL                  negative values of parcel water
!LL                  content are not generated.
!LL                        Pete Inness.
!LL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!     6.2  03/02/05  Added section 5A. R A Stratton
!     6.4  12/12/06  Removed def 3C.  R A Stratton
!
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
!LL  VERSION NO. 1
!LL
!LL  LOGICAL COMPONENTS COVERED : P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE THETAR (NPNTS,THRK,QRK,XSQR,BGMK,THEK,QEK,QPK,QSEK,    &
     &                   DQSK,BWKP1,EXK,PK)
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
!
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NPNTS          ! VECTOR LENGTH
!
      INTEGER I              ! LOOP COUNTER
!
!
!----------------------------------------------------------------------
! VARIABLES THAT ARE INPUT
!----------------------------------------------------------------------
!
      REAL THEK(NPNTS)       ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                             !    IN LAYER K (K)
!
      REAL QEK(NPNTS)        ! IN ENVIRONMENT MIXING RATIO
                             !    IN LAYER K (KG/KG)
!
      REAL QPK(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K
                             !    (KG/KG)
!
      REAL QSEK(NPNTS)       ! IN SATURATION MIXING RATIO OF THE
                             !    ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL DQSK(NPNTS)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    ENVIRONMENT OF LAYER K (KG/KG/K)
!
      LOGICAL BGMK(NPNTS)    ! IN MASK FOR PARCELS SATURATED IN LAYER K
!
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR POINTS AT WHICH CONDENSATE
                             !    IS LIQUID IN LAYER K+1
!
      REAL EXK(NPNTS)        ! IN EXNER RATIO FOR LEVEL K
!
      REAL PK(NPNTS)         ! IN PRESSURE AT LEVEL K (PA)
!
!
!----------------------------------------------------------------------
! VARIABLES THAT ARE OUTPUT
!----------------------------------------------------------------------
!
      REAL THRK(NPNTS)       ! OUT PARCEL DETRAINMENT POTENTIAL
                             !     TEMPERATURE IN LAYER K (K)
!
      REAL QRK(NPNTS)        ! OUT PARCEL DETRAINMENT MIXING RATIO
                             !     IN LAYER K (KG/KG)
!
      REAL XSQR(NPNTS)       ! OUT EXCESS WATER VAPOUR OF
                             !     DETRAINING AIR (KG/KG)
!
!
!----------------------------------------------------------------------
! VARIABLES THAT ARE DEFINED LOCALLY
!
      REAL TT(NPNTS)         ! TEMPORARY TEMPERATURE FOR CALCULATION
                             ! OF SATURATION MIXING RATIO (K)
!
!
!----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!----------------------------------------------------------------------
!
      EXTERNAL QSAT
!
!*----------------------------------------------------------------------
!
      DO 20 I=1,NPNTS
!L
!L----------------------------------------------------------------------
!L  CALCULATE THE POTENTIAL TEMPERATURE OF DETRAINING AIR
!L
!L  UM DOCUMENTATION PAPER P27
!L  SECTION (6), EQUATION (26)
!L----------------------------------------------------------------------
!L
       IF (.NOT.BGMK(I)) THEN
          THRK(I)=THEK(I) * (1. + C_VIRTUAL*QEK(I)) /                   &
     &                 (1. + C_VIRTUAL*QPK(I))
       ELSE
          THRK(I) = THEK(I)*(1.0 + C_VIRTUAL*(QEK(I)-QSEK(I))/          &
     &                 (1.0 + C_VIRTUAL*THEK(I)*DQSK(I)))
       ENDIF
!L
!L----------------------------------------------------------------------
!L  CALCULATE THE MIXING RATIO OF THE DETRAINING AIR AIR THE
!L  DIFFERENCE BETWEEN THIS AND THE MIXING RATIO OF THE MEAN
!L  PARCEL IN LAYER K
!L
!L  THE MOISTURE DIFFERENCE IS USED TO CALCULATE THE
!L  COND_DET_K TERM OF EQUATION (30), SECTION (6),
!L  UM DOCUMENTATIONM PAPER P27
!L----------------------------------------------------------------------
!L
!
!-----------------------------------------------------------------------
! CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
! PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
! MIXING RATIO
!-----------------------------------------------------------------------
!
       TT(I) = THRK(I)*EXK(I)
   20  CONTINUE
! DEPENDS ON: qsat
      CALL QSAT (XSQR,TT,PK,NPNTS)
!
      DO 30 I=1,NPNTS
!L----------------------------------------------------------------------
!L  SMALL NUMERICAL APPROXIMATIONS IN THE ABOVE CALCULATIONS CAN MEAN
!L  THAT THE DETRAINING PARCEL IS NO LONGER SATURATED AT THRK. ADD A
!L  CHECK TO SEE IF THE PARCEL IS STILL SATURATED, AND RESET BGMK TO
!L  FALSE IF IT IS NOT.
!L---------------------------------------------------------------------
         IF(XSQR(I) >  QPK(I))BGMK(I)=.FALSE.


       IF (BGMK(I)) THEN
          QRK(I)  = XSQR(I)
          XSQR(I) = QPK(I) - XSQR(I)
       ELSE
          QRK(I)  = QPK(I)
          XSQR(I) = 0.
       ENDIF
   30  CONTINUE
!
      RETURN
      END SUBROUTINE THETAR
