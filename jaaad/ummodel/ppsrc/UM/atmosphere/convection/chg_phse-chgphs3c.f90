
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE CHG_PHSE----------------------------------------------
!LL
!LL  PURPOSE : CHANGE OF PHASE ROUTINE FOR POINTS WHERE NO
!LL            DOWNDRAUGHT OCCURING
!LL
!LL            UPDATES POTENTIAL TEMPERATURE OF LAYER K
!LL            AS PRECIPITATION CHANGES PHASE IN SITU
!LL
!LL            ADD LATENT HEATING WHERE PRECIPITATION CROSSES A
!LL            MELTING OR FREEZING LEVEL
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL  VERSION NO. 4  DATED 5/2/92
!LL
!LL  PROJECT TASK : P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE CHG_PHSE (NPNTS,K,RAIN,SNOW,DTHBYDT_KM1,               &
                           EXK,EXKM1,DELPKM1,THE_K,THE_KM1,TIMESTEP,    &
                           CCA)

      Use cv_cntl_mod, Only:                                            &
          lcv_phase_lim

      IMPLICIT NONE
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
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
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
! CLDAREA start

      ! fractional cloud area when no dd, used in evaporation calc
      REAL,PARAMETER:: CLDAREA=1.0

! CLDAREA end
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NPNTS                ! IN VECTOR LENGTH
!
      INTEGER I                    ! LOOP COUNTER
!
      INTEGER K                    ! IN MODEL LAYER
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      REAL EXK(NPNTS)              ! IN EXNER RATIO FOR LAYER K
!
      REAL EXKM1(NPNTS)            ! IN EXNER RATIO FOR LAYER K-1
!
      REAL DELPKM1(NPNTS)          ! IN PRESSURE DIFFERENCE ACROSS
                                   !    LAYER K-1 (PA)
!
      REAL THE_K(NPNTS)            ! IN POTENTIAL TEMPERATURE OF
                                   !    ENVIRONMENT IN LAYER K
!
      REAL THE_KM1(NPNTS)          ! IN POTENTIAL TEMPERATURE OF
                                   !    ENVIRONMENT IN LAYER K-1
!
      REAL TIMESTEP                ! IN Timestep in seconds
!
      REAL CCA(NPNTS)              ! IN Convective cloud area

!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!----------------------------------------------------------------------
!
      REAL RAIN(NPNTS)             ! INOUT
                                   ! IN  AMOUNT OF FALLING RAIN
                                   !     (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     RAIN (KG/M**2/S)
!
      REAL SNOW(NPNTS)             ! INOUT
                                   ! IN  AMOUNT OF FALLING SNOW
                                   !     (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     SNOW (KG/M**2/S)
!
      REAL DTHBYDT_KM1(NPNTS)      ! INOUT
                                   ! IN  INCREMENT TO MODEL POTENTIAL
                                   !     TEMPERATURE IN LAYER K-1
                                   ! OUT UPDATED INCREMENT TO MODEL
                                   !     POTENTIAL TEMPERATURE IN LAYER
                                   !     K-1 DUE TO CHANGE OF PHASE
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!---------------------------------------------------------------------
!
      REAL FACTOR                  ! USED IN THE CALCULATION OF
                                   ! CHANGE OF PHASE OF FALLING
                                   ! PRECIPITATION
!
      REAL WPC                     ! AMOUNT OF PRECIPITATION WHICH CAN
!                                  ! CHANGE PHASE
!
      REAL CA                      ! LOCAL CLOUD AREA
!
      REAL THE_KM1_NEW             ! THE_KM1 UPDATED WITH ALL INCREMENTS
!                                  ! PRIOR TO THIS SUBROUTINE
!
      LOGICAL BPPNWT_K             ! MASK WHERE PRECIPITATION IS LIQUID
                                   ! IN LAYER K
!
      LOGICAL BPPNWT_KM1           ! MASK WHERE PRECIPITATION IS LIQUID
                                   ! IN LAYER K-1
!
!L
!L----------------------------------------------------------------------
!L  ADD LATENT HEATING WHERE PRECIP CROSSES A MELTING OR FREEZING LEVEL
!L
!L  UM DOCUMENTATION PAPER 27
!L  SECTION (11), EQUATION (42)
!L----------------------------------------------------------------------
!L
!
      IF (LCV_PHASE_LIM) THEN
!
      DO I=1,NPNTS
          THE_KM1_NEW = THE_KM1(I) + TIMESTEP*DTHBYDT_KM1(I)
          BPPNWT_K = THE_K(I) >  TM/EXK(I)
          BPPNWT_KM1 = THE_KM1_NEW  >   TM/EXKM1(I)
          FACTOR = LF*G/(EXKM1(I)*CP*DELPKM1(I))
          CA = CLDAREA * CCA(I)
!
! FREEZE
!
          IF (.NOT.BPPNWT_KM1.AND.(BPPNWT_K.OR.RAIN(I) >  0.0)) THEN
            WPC = MIN( RAIN(I),                                         &
     &               CA * (TM/EXKM1(I) - THE_KM1_NEW) /                 &
     &                    (TIMESTEP * FACTOR) )
            DTHBYDT_KM1(I) = DTHBYDT_KM1(I) + WPC * FACTOR
            SNOW(I) = SNOW(I) + WPC
            RAIN(I) = RAIN(I) - WPC
          END IF
!
! MELT
!
          IF (BPPNWT_KM1.AND.(.NOT.BPPNWT_K.OR.SNOW(I) >  0.0)) THEN
            WPC = MIN( SNOW(I),                                         &
     &               CA * (THE_KM1_NEW - TM/EXKM1(I)) /                 &
     &                    (TIMESTEP * FACTOR) )
            DTHBYDT_KM1(I) = DTHBYDT_KM1(I) - WPC * FACTOR
            RAIN(I) = RAIN(I) + WPC
            SNOW(I) = SNOW(I) - WPC
          END IF
        END DO
!
      ELSE
!
        DO I=1,NPNTS
        BPPNWT_K = THE_K(I)*EXK(I) >  TM
        BPPNWT_KM1 = THE_KM1(I)*EXKM1(I) >  TM
        FACTOR = LF*G/(EXKM1(I)*CP*DELPKM1(I))
! FREEZE
        IF (.NOT.BPPNWT_KM1.AND.(BPPNWT_K.OR.RAIN(I) >  0.0)) THEN
           DTHBYDT_KM1(I) = DTHBYDT_KM1(I)+RAIN(I)*FACTOR
           SNOW(I) = SNOW(I)+RAIN(I)
           RAIN(I) = 0.0
        END IF
! MELT
        IF (BPPNWT_KM1.AND.(.NOT.BPPNWT_K.OR.SNOW(I) >  0.0)) THEN
           DTHBYDT_KM1(I) = DTHBYDT_KM1(I)-SNOW(I)*FACTOR
           RAIN(I) = RAIN(I)+SNOW(I)
           SNOW(I) = 0.0
        END IF
      END DO
!
      ENDIF     ! IF LCV_PHASE_LIM
!
      RETURN
      END SUBROUTINE CHG_PHSE
!
