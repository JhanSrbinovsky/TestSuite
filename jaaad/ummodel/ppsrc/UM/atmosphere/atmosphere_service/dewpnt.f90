
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINE DEWPNT-------------------------------------------------
!LL
!LL  Purpose: Calculates the 1.5 metre dewpoint from 1.5 metre specific
!LL           humidity, 1.5 metre temperature and 1.5 metre pressure.
!LL
!LL  Suitable for single column usage.
!LL
!LL  Model            Modification history:
!LL version  Date
!LL
!LL    3.3  28/04/94 Created by Steve Woltering
!LL    4.4  Sept 97  Avoid crash if negative Q input. Damian Wilson.
!      5.5 17/04/03 Remove reference to obsolete section
!                   C90_1A. T.White
!      6.2  03/02/06 Move to a71_1a. P.Selwood
!LL
!LL  Programming standard:  Unified Model Documentation Paper No 3,
!LL                         Version 5, dated 08/12/92
!LL Documentation:  To be added to UM Doc Paper ?
!LL
!LLEND-----------------------------------------------------------------
!
!*L
!*LArguments:----------------------------------------------------------
      SUBROUTINE DEWPNT(                                                &
     & Q, P, T,                                                         &
                     ! IN
     & P_FIELD,                                                         &
                     ! IN
     & TD                                                               &
                     ! OUT
     &)
      IMPLICIT NONE
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

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
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
      INTEGER P_FIELD         ! IN Size of field arrays.
      REAL P(P_FIELD),                                                  &
                              ! IN Pressure.
     &     Q(P_FIELD),                                                  &
                              ! IN Specific humidity.
     &     T(P_FIELD)         ! IN Temperature.
      REAL RV,                                                          &
                              ! LOCAL Gas constant for water vapour.
     &     RL1,                                                         &
                              ! LOCAL Latent heat of evaporation.
     &     RT,                                                          &
                              ! LOCAL.
     &     P1(P_FIELD),                                                 &
                              ! LOCAL Pressure.
!                               j/Kg at 0 deg C.
     &     RL(P_FIELD),                                                 &
                              ! LOCAL.
     &     Q0(P_FIELD),                                                 &
                              ! LOCAL local SH.
     &     ES0,                                                         &
                              ! LOCAL Saturated vapour pressure.
     &     V_PRES(P_FIELD)    ! LOCAL Vapour pressure.
      INTEGER I               ! LOCAL loop variable.
      REAL TD(P_FIELD)        ! OUT Dew point.
      PARAMETER ( RV = R / EPSILON )
      PARAMETER ( RL1 = -2.73E3 )
!*----------------------------------------------------------------------
!*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
      EXTERNAL  QSAT_WAT
!----------------------------------------------------------------------
!  Calculate P in HPa.
!
      DO I=1,P_FIELD
        P1(I) = P(I) / 100.0
!----------------------------------------------------------------------
!  Calculate RL - The latent heat of evaporation.
        RL(I) = LC + RL1 * ( T(I) - TM )
!----------------------------------------------------------------------
!  Calculate Vapour pressure, and from that the dewpoint in Kelvins.
        V_PRES(I) = Q(I) * P1(I) / ( EPSILON + Q(I))
      ENDDO
! DEPENDS ON: qsat_wat
      CALL QSAT_WAT(Q0,T,P,P_FIELD)
      DO I=1,P_FIELD
        IF (V_PRES(I)  >   0.0) THEN
          ES0=(Q0(I) * P1(I)) / (EPSILON + Q0(I))
          RT = (1 / T(I)) - ( RV * ALOG(V_PRES(I)/ES0) )/RL(I)
          TD(I)=1.0/RT
          IF (TD(I)  >   T(I)) TD(I) = T(I)
        ELSE
          TD(I)=0.0
!         print*,'WARNING. Neg or zero Q in dewpoint calc.'
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DEWPNT
