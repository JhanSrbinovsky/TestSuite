
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
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
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
