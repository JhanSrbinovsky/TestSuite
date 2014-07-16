
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DET_RATE-----------------------------------------------
!LL
!LL  PURPOSE : CALCULATES THE FORCED DETRAINMENT RATE IN LAYER K
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL   4.4   17/10/97  New version optimised for T3E.
!LL                   Single PE optimisations           D.Salmond
!     6.2   03/02/05  Added section 5A R A Stratton
!     6.4   12/12/06  Removed def 3C. R A Stratton
!LL
!LL  PROGRAMMING STANDARDS :
!LL
!LL  LOGICAL COMPONENTS COVERED: P27
!LL
!LL  SYSTEM TASK :
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL                  SECTION (6), EQUATION (31)
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE DET_RATE (NPNTS,DELTAK,THRK,XSQR,THPK,THEK,THEKP1,     &
     &                   XSQKP1,THPKP1,BWKP1,BCALC,EKP14,EKP34,         &
     &                   EXK,EXKP1)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS
!-----------------------------------------------------------------------
!
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
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
      INTEGER NPNTS            ! VECTOR LENGTH
!
      INTEGER I                ! LOOP COUNTER
!
!
!-----------------------------------------------------------------------
! VARIABLES THAT ARE INPUT
!-----------------------------------------------------------------------
!
      REAL THRK(NPNTS)         ! IN PARCEL DETRAINMENT POTENTIAL
                               !    TEMPERATURE IN LAYER K (K)
!
      REAL XSQR(NPNTS)         ! IN EXCESS WATER VAPOUR OF THE
                               !    DETRAINING AIR IN LAYER K (KG/KG)
!
      REAL THPK(NPNTS)         ! IN PARCEL POTENTIAL TEMPERATURE
                               !    IN LAYER K (K)
!
      REAL THEK(NPNTS)         ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                               !    IN LAYER K (K)
!
      REAL THEKP1(NPNTS)       ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                               !    IN LAYER K+1 (K)
!
      REAL XSQKP1(NPNTS)       ! IN EXCESS WATER VAPOUR OF THE PARCEL
                               !    IN LAYER K+1 (KG/KG)
!
      REAL THPKP1(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE
                               !    IN LAYER K+1 (K)
!
      LOGICAL BCALC(NPNTS)     ! IN MASK FOR POINTS AT WHICH
                               !    CALCULATIONS OF THIS ROUTINE
                               !    ARE NEEDED
!
      LOGICAL BWKP1(NPNTS)     ! IN MASK FOR THOSE POINTS AT WHICH
                               !    CONDENSATE IS LIQUID IN LAYER K+1
!
      REAL EKP14(NPNTS)        ! IN ENTRAINEMNT RATE FOR LEVEL K+1/4
                               !    MULTIPLIED BY APPROPRIATE LAYER
                               !    THICKNESS
!
      REAL EKP34(NPNTS)        ! IN ENTRAINEMNT RATE FOR LEVEL K+3/4
                               !    MULTIPLIED BY APPROPRIATE LAYER
                               !    THICKNESS
!
      REAL EXK(NPNTS)          ! IN EXNER RATIO FOR LEVEL K
!
      REAL EXKP1(NPNTS)        ! IN EXNER RATIO FOR LEVEL K+1
!
!
!-----------------------------------------------------------------------
! VARIABLES THAT ARE OUTPUT
!-----------------------------------------------------------------------
!
      REAL DELTAK(NPNTS)       ! OUT PARCEL FORCED DETRAINMENT RATE
                               !     IN LAYER K MULTIPLIED BY
                               !     APPROPRIATE LAYER THICKNESS
!
!
!-----------------------------------------------------------------------
! VARIABLES THAT ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
      REAL EL                  ! LATENT HEAT OF CONDENSATION OR
                               ! (CONDENSATION + FUSION) (J/KG)
!
      REAL EPSS                ! (1+EKP14)*(1+EKP34)
!
!*---------------------------------------------------------------------
!L
!L---------------------------------------------------------------------
!L  NO SIGNIFICANT STRUCTURE
!L---------------------------------------------------------------------
!L
!
      DO 10 I=1,NPNTS
       IF (BCALC(I)) THEN
!
!-----------------------------------------------------------------------
!   CREATE A VECTOR OF LATENT HEATS
!-----------------------------------------------------------------------
!
       IF (BWKP1(I)) THEN
          EL = LC
       ELSE
          EL = LC + LF
       ENDIF
!
!-----------------------------------------------------------------------
!   CALCULATE DETRAINMENT RATES
!-----------------------------------------------------------------------
!
       EPSS = (1. + EKP14(I)) * (1. + EKP34(I))
          DELTAK(I) = EKP14(I)*THEK(I)                                  &
     &        + EKP34(I)*(1.+EKP14(I))*THEKP1(I)                        &
     &        - EPSS*(THPKP1(I) - EL/(EXKP1(I)*CP) * XSQKP1(I))
!
          DELTAK(I) =   (DELTAK(I) + THPK(I))                           &
     &  *(EXK(I)*CP)/((DELTAK(I) + THRK(I))*(EXK(I)*CP) - EL*XSQR(I))
!
!----------------------------------------------------------------------
!  FROM A THEORETICAL VIEW POINT DELTAK CANNOT = 1 . HOWEVER
!  BECAUSE OF APPROXIMATION USED IN THE CALCULATION NUMERICALLY IT
!  MAY BE POSSIBLE.  HENCE IF DELTAK = 1 SET IT TO SLIGHTLY SMALLER
!  THAN 1
!----------------------------------------------------------------------
!
          DELTAK(I) = MIN(0.99999,DELTAK(I))
!
       ENDIF
   10 CONTINUE
!
      RETURN
      END SUBROUTINE DET_RATE
