
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE BOUY_TQ

! PURPOSE: To calculate buoyancy parameters on p,T,q-levels
!
! METHOD:
!
! HISTORY:
! DATE   VERSION   COMMENT
! ----   -------   -------
! new deck
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!LL   5.2  27/09/00   change from QSAT_2D to QSAT           A.Malcolm
!!!  5.3    Feb. 01  Calculate grid-box mean parameters (APLock)
!  6.1  17/05/04  Change Q, QCL, QCF dims to enable substepping.
!                                                       M. Diamantakis
!
! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
!   THIS CODE IS WRITTEN TO UMDP 3 PROGRAMMING STANDARDS.
!

      SUBROUTINE BOUY_TQ (                                              &
     & row_length, rows, halo_i,halo_j,BL_LEVELS, LQ_MIX_BL             &
     &,P,T,Q,QCF,QCL,CF                                                 &
     &,BT,BQ,BT_CLD,BQ_CLD,BT_GB,BQ_GB,A_QS,A_DQSDT,DQSDT               &
     &,LTIMER                                                           &
     & )

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL LTIMER          ! IN Flag for TIMER diagnostics


      LOGICAL LQ_MIX_BL       ! IN switch for using mixing ratios
      INTEGER                                                           &
     & row_length, rows, halo_i,halo_j                                  &
     &,BL_LEVELS              ! IN No. of atmospheric levels for which
!                                boundary layer fluxes are calculated.
!                                Assumed  <=30 for dimensioning GAMMA()
!                                in common deck C_GAMMA

      REAL                                                              &
     & P(row_length,rows,BL_LEVELS)                                     &
                                      ! IN Pressure at pressure points.
     &,T(row_length,rows,BL_LEVELS)                                     &
                                      ! IN Temperature (K). At P points
     &,Q(row_length,rows,BL_LEVELS)                                     &
                                ! IN Sp humidity (kg water per kg air).
     &,QCL(row_length,rows,BL_LEVELS)                                   &
                                ! IN Cloud liq water (kg per kg air).
     &,QCF(row_length,rows,BL_LEVELS)                                   &
                                ! IN Cloud liq water (kg per kg air).
     &,CF(row_length, rows, BL_LEVELS)! IN Cloud fraction (decimal).


! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

      REAL                                                              &
     & BQ(row_length,rows,BL_LEVELS)                                    &
                                ! OUT A buoyancy parameter for clear air
     &,BT(row_length,rows,BL_LEVELS)                                    &
                                ! OUT A buoyancy parameter for clear air
     &,BQ_CLD(row_length,rows,BL_LEVELS)                                &
!                             ! OUT A buoyancy parameter for cloudy air
     &,BT_CLD(row_length,rows,BL_LEVELS)                                &
!                             ! OUT A buoyancy parameter for cloudy air
     &,BQ_GB(row_length,rows,BL_LEVELS)                                 &
                                ! OUT A grid-box mean buoyancy parameter
     &,BT_GB(row_length,rows,BL_LEVELS)                                 &
                                ! OUT A grid-box mean buoyancy parameter
     &,A_QS(row_length,rows,BL_LEVELS)                                  &
!                             ! OUT Saturated lapse rate factor
     &,A_DQSDT(row_length,rows,BL_LEVELS)                               &
!                             ! OUT Saturated lapse rate factor
     &,DQSDT(row_length,rows,BL_LEVELS)
!                             ! OUT Derivative of q_SAT w.r.t. T

! LOCAL VARIABLES.

      REAL                                                              &
     & QS(row_length,rows)            ! WORK Saturated mixing ratio.

      INTEGER                                                           &
     &  I,j                                                             &
     &, K

      REAL                                                              &
     &  BC

      EXTERNAL                                                          &
     &  QSAT, TIMER

!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
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
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
! C_SOILH start
      ! No. of soil layers (must = NSOIL).
      REAL,PARAMETER:: PSOIL=4

      ! Tunable characteristic freq (rad/s)
      REAL,PARAMETER:: OMEGA1=3.55088E-4

      ! Density of lying snow (kg per m**3)
      REAL,PARAMETER:: RHO_SNOW=250.0

      ! Depth of `effective' snow surface layer (m)
      REAL,PARAMETER:: DEFF_SNOW=0.1

      ! Thermal conductivity of lying snow (Watts per m per K).
      REAL,PARAMETER:: SNOW_HCON=0.265

      ! Thermal capacity of lying snow (J/K/m3)
      REAL,PARAMETER:: SNOW_HCAP=0.63E6

! C_SOILH end
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------


      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP
      PARAMETER (                                                       &
     & ETAR=1.0/(1.0-EPSILON)                                           &
                                ! Used in buoyancy parameter BETAC.
     &,GRCP=G/CP                                                        &
                                ! Used in DZTL, FTL calculations.
     &,LCRCP=LC/CP                                                      &
                                ! Latent heat of condensation / CP.
     &,LFRCP=LF/CP                                                      &
                                ! Latent heat of fusion / CP.
     &,LS=LC+LF                                                         &
                                ! Latent heat of sublimation.
     &,LSRCP=LS/CP                                                      &
                                ! Latent heat of sublimation / CP.
     &)

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BOUY_TQ ',3)
      ENDIF
!-----------------------------------------------------------------------
!! 1.  Loop round levels.
!-----------------------------------------------------------------------
      DO K=1,BL_LEVELS
!-----------------------------------------------------------------------
!! 1.1 Calculate saturated specific humidity at pressure and
!!     temperature of current level.
!-----------------------------------------------------------------------
! DEPENDS ON: qsat_mix
        CALL QSAT_mix(QS,T(1,1,K),P(1,1,K),row_length*rows,Lq_mix_bl)
!
        do j=1,rows
        DO I=1,row_length
!ajm        DO I=P1,P1+P_POINTS-1

!-----------------------------------------------------------------------
!! 1.2 Calculate buoyancy parameters BT and BQ, required for the
!!     calculation of stability.
!-----------------------------------------------------------------------

          BT(I,j,K) = 1.0/T(I,j,K)
          BQ(I,j,K) =                                                   &
     &      C_VIRTUAL/(1.0+C_VIRTUAL*Q(I,j,K)-QCL(I,j,K)-QCF(I,j,K))
!
          IF (T(I,j,K)  >   TM) THEN
            DQSDT(I,j,K) = (EPSILON * LC * QS(I,j))                     &
     &                   / ( R * T(I,j,K) * T(I,j,K) )
!                      ...  (Clausius-Clapeyron) for T above freezing
!
            A_QS(I,j,K) = 1.0 / (1.0 + LCRCP*DQSDT(I,j,K))
!
            A_DQSDT(I,j,K) = A_QS(I,j,K) * DQSDT(I,j,K)
!
            BC = LCRCP*BT(I,j,K) - ETAR*BQ(I,j,K)
!
          ELSE
            DQSDT(I,j,K) = (EPSILON * LS * QS(I,j))                     &
     &                   / ( R * T(I,j,K) * T(I,j,K) )
!                      ...  (Clausius-Clapeyron) for T below freezing
!
            A_QS(I,j,K) = 1.0 / (1.0 + LSRCP*DQSDT(I,j,K))
!
            A_DQSDT(I,j,K) = A_QS(I,j,K) * DQSDT(I,j,K)
!
            BC = LSRCP*BT(I,j,K) - ETAR*BQ(I,j,K)
!
          ENDIF
!
!-----------------------------------------------------------------------
!! 1.3 Calculate in-cloud buoyancy parameters.
!-----------------------------------------------------------------------
!
          BT_CLD(I,j,K) = BT(I,j,K) - A_DQSDT(I,j,K) * BC
          BQ_CLD(I,j,K) = BQ(I,j,K) + A_QS(I,j,K) * BC

!-----------------------------------------------------------------------
!! 1.4 Calculate grid-box mean buoyancy parameters.
!-----------------------------------------------------------------------
!
          BT_GB(I,j,K) = BT(I,j,K) +                                    &
     &                   CF(I,j,K)*( BT_CLD(I,j,K) - BT(I,j,K) )
          BQ_GB(I,j,K) = BQ(I,j,K) +                                    &
     &                   CF(I,j,K)*( BQ_CLD(I,j,K) - BQ(I,j,K) )
!
        ENDDO ! p_points,j
        ENDDO ! p_points,i
      ENDDO ! bl_levels

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BOUY_TQ ',4)
      ENDIF
      RETURN
      END SUBROUTINE BOUY_TQ
