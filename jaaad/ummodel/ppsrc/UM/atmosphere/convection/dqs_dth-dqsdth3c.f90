
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DQSDTH-------------------------------------------------
!LL
!LL  PURPOSE : CALCULATES GARDIENT OF SATURATION MIXING RATIO
!LL            WITH POTENTIAL TEMPERATURE FORM THE
!LL            CLAUSIUS-CLAPEYRON EQUATION
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL   4.4   17/10/97  New version optimised for T3E.
!LL                   Removed divide
!LL                                                   D.Salmond
!     6.2   03/02/05  Added section 5A. R A Stratton
!     6.4   12/12/06  Removed def 3C,. R A Stratton
!LL
!LL  PROGRAMMING STANDARDS :
!LL
!LL  LOGICAL COMPONENTS COVERED: P27
!LL
!LL  SYSTEM TASK :
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL                  SECTION(4), EQUATION (20)
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE DQS_DTH (DQS,K,THEK,QSEK,EXK,NPNTS)
      IMPLICIT NONE
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
!
! RV gas constant for water vapour (J/kg/K)
      REAL,PARAMETER:: RV = 461.1
! RV end
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NPNTS        ! IN VECTOR LENGTH
!
      INTEGER K            ! IN PRESENT MODEL LAYER
!
      INTEGER I            ! LOOP COUNTER
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      REAL THEK(NPNTS)     ! IN POTENTIAL TEMPERATURE FOR LAYER K (K)
!
      REAL QSEK(NPNTS)     ! IN SATURATION MIXING RATIO FOR LAYER K (K)
!
      REAL EXK(NPNTS)      ! IN EXNER RATIO FOR LAYER K
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!----------------------------------------------------------------------
!
      REAL DQS(NPNTS)      ! OUT GRADIENT OF SATURATION MIXING RATIO
                           !     WITH POTENTIAL TEMPERATURE FOR LAYER K
                           !     (KG/KG/S)
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
!
      REAL EL              ! LATENT HEATING OF CONDENSATION OR
                           ! (CONDENSATION PLUS HEATING) (J/KG)
!
!*---------------------------------------------------------------------
!L
!L---------------------------------------------------------------------
!L NO SIGNIFICANT STRUCTURE
!L---------------------------------------------------------------------
!L
      DO 10 I=1,NPNTS
!
!-----------------------------------------------------------------------
!  CREATE A VECTOR OF LATENT HEATS ACCORDING TO WHETHER QSAT IS WRT
!  ICE OR WATER
!-----------------------------------------------------------------------
!
      IF (THEK(I)*EXK(I)  >   273.) THEN
          EL = LC
       ELSE
          EL = LC + LF
       ENDIF
!
!-----------------------------------------------------------------------
! CALCULATE D(QSAT)/D(THETA)
!-----------------------------------------------------------------------
!
       DQS(I) = EL*QSEK(I)/(EXK(I)*RV*THEK(I)*THEK(I))
   10  CONTINUE
!
      RETURN
      END SUBROUTINE DQS_DTH
