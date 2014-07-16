
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE CON_RAD------------------------------------------------
!LL
!LL  PURPOSE : CALCULATES CONVECTIVE CLOUD TOP, BASE AND
!LL            AMOUNT
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL  3.3  23/12/93 Change to cloud top because of
!LL                change to detrainment rate calculation.  D.Gregory.
!LL
!LL  3.4  21/03/94  Add lowest conv.cloud diagnostics.  R.T.H.Barnes.
!LL
!LL  4.4  26/09/97  Pass in extra cloud water variable to allow rain
!LL                 out in CLOUDW before calculation of water path
!LL                 if L_CCW is set to .TRUE. in CLOUDW.      J.M.G.
!
!    5.4  20/08/02: Changes to logic to allow for a parcel starting to
!                   convect with non-zero (environmental value) cloud
!                   condensate.
!                                                    A.C.Bushell.
!    5.5  17/04/03 Removal of references to obsolete sections
!                  A05_2A,2C,3B. T.White
!
!    6.1  27/05/04: Change of A05_3C only compiler directive to PC2
!                   if test for 4A convection.   A. Kerr-Munslow
!    6.2  03/02/05  Added section 5A R.A.Stratton
!    6.4  12/12/06  Removed Def 3C. R A Stratton
!
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
!LL  VERSION NO. 1
!LL
!LL  LOGICAL COMPONENT NUMBER: P27
!LL
!LL  SYSTEM TASK : P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
!LL                  SECTION (9)
!LL
!LLEND-----------------------------------------------------------------
!
!
!+ Convective Cloud Top, Base and Amount Calculation Scheme.
! Subroutine Interface:
      SUBROUTINE CON_RAD (                                              &
     &           K,XPK,XPKP1,FLXKP1,BTERM,CCA,ICCB,ICCT,START_LEV,      &
     &           TCW,CCW,CCLWP,DELPKP1,LCCA,LCBASE,LCTOP,LCCLWP,NPNTS,  &
     &           L_Q_INTERACT)
!
      IMPLICIT NONE
!
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
!
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!
!----------------------------------------------------------------------
! VECTOR LENGTH AND LOOP VARIABLES
!----------------------------------------------------------------------
!
      INTEGER NPNTS        ! IN VECTOR LENGTH
!
      INTEGER K            ! IN PRESENT MODEL LAYER
!
      INTEGER I            ! LOOP COUNTER
!
! ----------------------------------------------------------------------
! Arguments with Intent IN.  ie:  Input variables.
! ----------------------------------------------------------------------
!
      REAL XPK(NPNTS)      ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
!
      REAL XPKP1(NPNTS)    ! IN PARCEL CLOUD WATER IN LAYER K+1 (KG/KG)
!
      LOGICAL BTERM(NPNTS) ! IN MASK FOR POINTS WHERE CONVECTION
                           !    IS ENDING
!
      REAL FLXKP1(NPNTS)   ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
!
      REAL DELPKP1(NPNTS)  ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
!
      REAL CCW(NPNTS)      ! IN PARCEL CLOUD WATER AS CALCULATED BEFORE
                           !    PRECIPITATION. LAYER K+1 (KG/KG)
!
!
      INTEGER START_LEV(NPNTS)  ! IN LEVEL AT WHICH CONVECTION INITIATED
!
      LOGICAL L_Q_INTERACT ! IN Switch (PC2) allows overwriting of
                           ! parcel variables
!
! ----------------------------------------------------------------------
! Arguments with Intent IN/OUT. ie: input variables changed on output.
! ----------------------------------------------------------------------
!
      REAL TCW(NPNTS)      ! INOUT
                           ! IN  TOTAL CONDENSED WATER SUMMED TO
                           !     LAYER K (KG/M**2/S)
                           ! OUT TOTAL CONDENSED WATER SUMMED TO
                           !     LAYER K+1 OR IF CONVECTION HAS
                           !     TERMINATED ZEROED (KG/M**2/S)
!
      REAL CCLWP(NPNTS)    ! INOUT
                           ! IN  TOTAL CLOUD LIQUID WATER PATH
                           !     SUMMED TO LAYER K  (KG/M**2)
                           ! OUT TOTAL CLOUD LIQUID WATER PATH
                           !     SUMMED TO LAYER K+1 (KG/M**2)
      REAL LCCA(NPNTS)      ! INOUT LOWEST CONV.CLOUD AMOUNT (%)
!
      INTEGER LCBASE(NPNTS) ! INOUT LOWEST CONV.CLOUD BASE LEVEL
!
      INTEGER LCTOP(NPNTS)  ! INOUT LOWEST CONV.CLOUD TOP LEVEL
!
! ----------------------------------------------------------------------
! Arguments with Intent OUT. ie: Output variables.
! ----------------------------------------------------------------------
!
      REAL CCA(NPNTS)      ! OUT CONVECTIVE CLOUD AMOUNT (%)
!
      INTEGER ICCB(NPNTS)   ! OUT CONVECTIVE CLOUD BASE LEVEL
!
      INTEGER ICCT(NPNTS)   ! OUT CONVECTIVE CLOUD TOP LEVEL
!
      REAL LCCLWP(NPNTS)    ! OUT LOWEST CONV.CLOUD LIQ.WATER PATH
!
!  External subroutine calls: ------------------------------------------
!
!     EXTERNAL None.
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
!
! ----------------------------------------------------------------------
!  CALCULATE Cloud Base and Lowest Cloud Base
!
!  WHEN CLOUD BASE SET ZERO TOTAL CONDENSED WATER
! ----------------------------------------------------------------------
!
      DO I=1, NPNTS
        IF ( XPK(I)  <=  0.0 .AND. CCW(I)  >   0.0 ) THEN
!       Assuming initial parcel condensate is zero
          ICCB(I)=K+1
          CCLWP(I)=0.0
!
          IF ( LCBASE(I)  ==  0 ) THEN
!         Lowest cloud base
            LCBASE(I) = K+1
            LCCLWP(I) = 0.0
          END IF  ! If_lcbase_1
!
        ELSE IF ( L_Q_INTERACT .AND. XPK(I)  >   0.0 .AND.              &
     &           K  ==  START_LEV(I) ) THEN
!       Non-zero initial parcel condensate (initialized to environment)
          ICCB(I)  = K
          CCLWP(I) = 0.0
!
          IF ( LCBASE(I)  ==  0 ) THEN
!         Lowest cloud base
            LCBASE(I) = K
            LCCLWP(I) = 0.0
          END IF  ! If_lcbase_2
!


        END IF
!
! ----------------------------------------------------------------------
!  CALCULATE Cloud Top and Lowest Cloud Top
! ----------------------------------------------------------------------
!
        IF (BTERM(I) .AND.                                              &
     &      ((CCW(I) >  0.0).OR.(XPK(I) >  0.0)) ) ICCT(I) = K+1

        IF (BTERM(I) .AND.  LCTOP(I) == 0 .AND.                         &
     &      ((CCW(I) >  0.0).OR.(XPK(I) >  0.0)) ) THEN
          LCTOP(I) = K+1
        END IF
!
        IF ( FLXKP1(I)  >   0.0) THEN
!L
!L---------------------------------------------------------------------
!L SUM TOTAL CONDENSED WATER PER SECOND - ASSUMES THAT THE INITIAL
!L CONVECTIVE LAYER IS UNSATURATED
!L---------------------------------------------------------------------
!L
          TCW(I) = TCW(I) + FLXKP1(I) * CCW(I) / G
!L
!L---------------------------------------------------------------------
!L SUM CONV CONDENSED WATER PATH - ASSUMES THAT THE INITIAL
!L CONVECTIVE LAYER IS UNSATURATED
!L---------------------------------------------------------------------
!L
          CCLWP(I) = CCLWP(I) + XPKP1(I) * DELPKP1(I) / G
!L
!L---------------------------------------------------------------------
!L SUM CONV CONDENSED WATER PATH up to lowest conv.cloud
!L ASSUMES THAT THE INITIAL CONVECTIVE LAYER IS UNSATURATED
!L---------------------------------------------------------------------
!L
          IF (LCCA(I) <= 0.0) THEN
            LCCLWP(I) = LCCLWP(I) + CCW(I) * DELPKP1(I) / G
          END IF
!
        END IF
!L
!L---------------------------------------------------------------------
!L CALCULATE CONVECTIVE CLOUD AMOUNT IF CONVECTION TERMINATES IN
!L LAYER K AND TOTAL CONDENSED WATER PATH OVER A TIME STEP
!L
!L UM DOCUMENTATION PAPER P27
!L SECTION (9), EQUATION (37)
!L---------------------------------------------------------------------
!L
        IF( BTERM(I) .AND. TCW(I) >  0.0 ) THEN
!
          IF ( TCW(I)  <   2.002E-6 ) TCW(I) = 2.002E-6
!
          CCA(I) = 0.7873 + 0.06 * LOG(TCW(I))
          IF (CCA(I)  >   1.0) CCA(I) = 1.0
!
          IF (LCCA(I) <= 0.0) THEN
            LCCA(I) = 0.7873 + 0.06 * LOG(TCW(I))
            IF (LCCA(I)  >   1.0) LCCA(I) = 1.0
          END IF
!
          TCW(I) = 0.0
!
        END IF
      END DO ! I loop over NPNTS
!
      RETURN
      END SUBROUTINE CON_RAD
