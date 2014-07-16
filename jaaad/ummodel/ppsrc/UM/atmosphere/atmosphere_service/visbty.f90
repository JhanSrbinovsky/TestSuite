
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE VISBTY -----------------------------------------------
!LL
!LL     PURPOSE:
!LL Process fields of temperature, specific humidity, cloud liquid
!LL water or ThetaL and qt to give visibility in metres.
!LL Calculated at model level (eg bottom eta level 25m)
!LL or level within surface layer eg screen ht ( 1.5M )
!LL
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.1  23/10/92  New deck author P.Smith
!LL  3.1  20/01/93  New deck - used as mods at 2.7 & 2.8/3.0.
!LL                 Interfacing done by R.T.H.Barnes.
!LL  3.2  29/04/93  CCN and derived constants moved to MODECK C_VISBTY
!LL                 Programmer Pete Clark.
!LL  3.4  07/06/94  Aerosol field introduced. Programmer Pete Clark.
!LL  4.0 05/09/95  Lower limit to aerosol introduced.  Pete Clark.
!    4.4  09/01/97  Only liquid water and not cloud ice is now used
!                   to calculate visibility. Damian Wilson.
!LL  4.5  29/04/98  Scheme replaced with NIMROD visibility diagnostic
!LL  5.0  22/11/99  5.0 introduced : pstar, ak, bk, replaced by
!LL                 p_layer.    JC Thil
!LL  5.0  14/02/00  Correction to prevent overwriting of D1 array
!LL                 when aerosol field updated with no aerosol in
!LL                 the model. R Rawlins
!    5.3  27/07/01  New code to take account of precip when calculating
!                   visibility.                              Pete Clark
!    5.5  17/04/03  Remove reference to obsolete section
!                   C90_1A. T.White
!    6.2  03/02/06  Moved to a71_1a. P.Selwood
!LL
!LL  Programming standard: U M Doc. Paper No. 4
!LL
!LL  Logical components covered :
!LL
!LL  Project task:
!LL
!LL  External documentation
!LL    Forecasting Research Scientific Paper NO.4
!LL    Diagnosis of visibility in the UK Met Office Mesoscale Model
!LL    and the use of a visibility analysis to constrain initial
!LL    conditions.  SP Ballard, BJ Wright, BW Golding    1992
!      NIMROD diagnostic:
!      Wright, B. J., 1997: Improvements to the Nimrod Visibility
!         Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!      Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!         for Nimrod. Met. Office FR Tech Rep., No. 222.
!LL
!LLEND----------------------------------------------------------------
!
!*L  Arguments:-------------------------------------------------------
      SUBROUTINE VISBTY(                                                &
     &             p_layer,T,Q,QCL,QCF                                  &
                                              !INPUT
     &           ,AEROSOL, PROB, RHCRIT, L_MURK                         &
                                                      !INPUT
     &           ,P_FIELD                                               &
                                                      !INPUT
     &           ,VISIBILITY)                         !OUTPUT
      IMPLICIT NONE
!---------------------------------------------------------------------
! Workspace usage:----------------------------------------------------
! 3 real arrays of size P_FIELD
!*--------------------------------------------------------------------
!*L-------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
      INTEGER                                                           &
     &        P_FIELD                   ! IN NO. points in field.
      REAL                                                              &
     &  p_layer(p_field)                                                &
     &       ,T(P_FIELD)                                                &
                                        ! IN Temperature
     &       ,Q(P_FIELD)                                                &
                                        ! IN Qt
     &       ,QCL(P_FIELD)                                              &
                                        ! IN cloud water array.
     &       ,QCF(P_FIELD)                                              &
                                        ! IN cloud ice array.
     &       ,AEROSOL(P_FIELD)                                          &
                                        ! IN Aerosol mixing ratio(ug/kg)
     &       ,PROB                                                      &
                                      ! IN Probability level ( e.g 0.5
                                      !    corresponds to median).
     &       ,RHCRIT                  ! IN Critical RH (determines
                                      !    width of distribiution)
      LOGICAL                                                           &
     &   L_MURK                        ! IN : Aerosol present
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
      REAL                                                              &
     &      VISIBILITY(P_FIELD)         ! OUT visibility array.
!*--------------------------------------------------------------------
!*L-------------------------------------------------------------------
! Local varables:-----------------------------------------------------
!---------------------------------------------------------------------
      REAL                                                              &
     &       QT(P_FIELD)                                                &
                                      ! total of cloud water and vapour
     &      ,P(P_FIELD)                                                 &
                                      ! pressure of level
     &      ,Qs(P_FIELD)              ! saturation vapour pressure
!*L  External subroutine called ----------------------------------------
      EXTERNAL QSAT_WAT
!*--------------------------------------------------------------------
! constants for visibility calculation used to be set here but now
! set in MODECK.
! PI needed to set new constants.
!---------------------------------------------------------------------
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
! C_DENSTY for subroutine SF_EXCH
      REAL,PARAMETER:: RHOSEA = 1026.0 ! density of sea water (kg/m3)
      REAL,PARAMETER:: RHO_WATER = 1000.0! density of pure water (kg/m3)
! C_DENSTY end
! C_VISBTY start
!LL Description:
!LL   This COMDECK contains declarations for constants used to diagnose
!LL visibility. Constants are set as PARAMTERs.
!LL
!LL
!LL  Model            Modification history:
!LL version  Date
!LL  3.2    29/04/93  CCN Parameters moved here from VISBTY so that
!LL                   they can also be used to compute fog fraction.
!LL                   Programmer: Pete Clark.
!LL  4.0 05/09/95  Variable AEROMAX used as upper limit to aerosol in
!LL                assimilation introduced. Programmer Pete Clark.
!LL  4.5 01/05/98  Completely re-written for NIMROD style diagnostic.
!    5.3 27/07/01  New parameters for the effect of precipitation
!                  visibility created.                 Pete Clark
!LL  5.3 17/10/01  Rename rho and rho_a. Adam Clayton
!LL
!LLEND----------------------------------------------------------------
      INTEGER, PARAMETER :: n_vis_thresh = 2

      ! Standard number density of the aerosol (/m3)
      REAL, PARAMETER :: N0 = 500.0E6

      ! Activation parameter
      REAL, PARAMETER :: B0= 0.5

      ! Radius of standard aerosol particle (m)
      REAL, PARAMETER :: radius0 = 0.16E-6

      REAL, PARAMETER :: FourThirds = 4.0/3.0     ! 4/3

      REAL :: vis_thresh(n_vis_thresh)
      DATA vis_thresh /1000.0,5000.0/
      ! Density of the the aerosol (Kg/m3)
      REAL, PARAMETER :: rho_aerosol = 1700.0

      ! Density of air (Kg/m3)
      REAL, PARAMETER :: rho_air = 1.0
      ! Standard aerosol mass mixing ratio (Kg/Kg)
      REAL, PARAMETER :: m0 = FourThirds * Pi *                         &
     &  radius0 * radius0 * radius0 *                                   &
     &                         (rho_aerosol/rho_air) * N0

      ! Aerosol particle radius/mass loading power
      REAL, PARAMETER :: power = 1.0/6.0

      ! Scattering coefficient normalisation
      REAL, PARAMETER :: Beta0  = 1.5 * Pi

      REAL, PARAMETER :: LiminalContrast   = 0.02

      ! Natural log of Liminal contrast
      REAL, PARAMETER :: LnLiminalContrast = -3.912023005

      ! Constant incorporating the scattering coefficient, normalisation
      ! transformation to visibility ( = ln(liminal contrast) / Beta0 )
      REAL, PARAMETER :: VisFactor = -LnLiminalContrast / Beta0

      ! Reciprocal of the clean air visibility
      REAL, PARAMETER :: RecipVisAir = 1.0E-5

      ! Constant involving surface energy of water
      REAL, PARAMETER :: A0 = 1.2E-9

      ! Visibility defining fog
      REAL, PARAMETER :: VISFOG = 1000.0

      ! Visibility defining mist
      REAL, PARAMETER :: VISMIST = 5000.0

      ! Minimum allowed aerosol
      REAL, PARAMETER :: AERO0 = 0.1

      ! maximum allowed aerosol
      REAL, PARAMETER :: AEROMAX = 200.0

      ! tunable parameter: the cumulative prob value at which vis is
      ! estimated
      REAL, PARAMETER :: calc_prob_of_vis = 0.4
!
! Local parameter variables
!
      REAL                                                              &
     &  OneThird                                                        &
                        ! 1/3
     &, RHmax                                                           &
                        ! Maximum value of relative humidity
                        !  which is allowed to feed into the
                        !  calculation of the 'fog' droplet radius
     &, RHmin                                                           &
                        ! Minimum value of relative humidity which
                        !  is allowed to feed into the calculation
                        !   of the 'fog' droplet radius
     &, Weight                                                          &
                        ! Weighting on new value for iterative
                        !  solution of droplet radius
     &, Delta_radius_star                                               &
                         ! Convergence required for iterative
                        !  solution of droplet radius
     &, N                                                               &
                        ! Local number density
     &, qt_limit                                                        &
                        ! Smallest Qt value allowed
     &, radius_star_min                                                 &
                        !
     &, radius_star_max                                                 &
                        !
     &, radius_star_factor!

      INTEGER                                                           &
     &  Niterations     !  Maximum number of iteration used to
                        !   estimate the water droplet radius
      PARAMETER ( OneThird = 1.0/3.0                                    &
     &,              RHmin = 0.001                                      &
     &,              RHmax = 0.99                                       &
     &,             Weight = 0.75                                       &
     &,  Delta_radius_star = 0.001                                      &
     &,        Niterations = 20                                         &
     &,           qt_limit = 0.0001                                     &
     &,    radius_star_min = 1.0                                        &
     &,    radius_star_max = 1000.0                                     &
     &, radius_star_factor = 4.0 )
!
! Local workspace variables
!
       INTEGER                                                          &
     &  Point                                                           &
                         !  Loop variable for points
     &, Iteration        !  Loop variable iterations used to estimate
                         !   the water droplet radius

      REAL                                                              &
     &  m_over_m0                                                       &
                         !  Ratio of  aerosol mass mixing ratio and
                         !   the standard aerosol mass mixing ratio
     &, RecipVis                                                        &
                         !  Recipirical of the visibility
     &, radius_dry                                                      &
                         !  Radius of dry aerosol particle (m)
     &, radius                                                          &
                         !  Radius of fog droplets (m)
     &, radius_star1                                                    &
                         !  Previous estimate of water droplet radius
                         !   divided by the dry radius
     &, radius_star2                                                    &
                         !  Current best estimate of water droplet
                         !   radius divided by the dry radius
     &, radius_act                                                      &
                         !  Activation droplet radius
     &, radius_star_act                                                 &
                         !  Activation droplet rad divided by dry rad
     &, A                                                               &
                         !  A0 divided by the dry radius
     &, RH_lim                                                          &
                         !  Limited RH value (fractional)
     &, Fn                                                              &
                         !  Value of droplet radius function
     &, Deriv                                                           &
                         !  Derivative of droplet radius function
     &, radius_star_diff                                                &
                         !  Absolute value of radius_star1 minus
                         !    radius_star2
     &, RHterm                                                          &
                         !  Relative humidity term in function to be
                         !   minimised to find the droplet radius
     &, qLterm                                                          &
                         !  Liquid water term in function to be
                         !   minimised to find the droplet radius
     &, RHderiv                                                         &
                         !  Derivative of relative humidity term
     &, qLderiv                                                         &
                         !  Derivative of liquid water term
     &, bs                                                              &
                         !  Width of distribution in total water
                         !   mixing ratio space (kg/kg)
     &, qt_mod                                                          &
                         !  Modified total water value based on the
                         !   probability of the value occurring
                         !   assuming a triangular distriubtion
                         !   of width bs.
     &, qt_mod_factor    !  Factor to multiply bs to modify qt
!     Check Prob is legal
      IF ( Prob  <   0.0 .OR.  Prob  >   1.0 ) THEN
        Write(6,*)"INVALID PROBABILITY VALUE in VISBTY",Prob
        Prob=MIN(MAX(Prob,0.0),1.0)
      ENDIF
!     Create factor to multiply bs by to modify qt
      IF ( Prob  ==  0.5 ) THEN
        qt_mod_factor = 0.0
      ELSE IF ( Prob  >=  0.0 .AND. Prob  <   0.5 ) THEN
        qt_mod_factor = ( 1.0 - SQRT( 2.0 * Prob ) )
      ELSE IF ( Prob  >=  0.5 .AND. Prob  <=  1.0 ) THEN
        qt_mod_factor = - ( 1.0 - SQRT( 2.0 * (1.0-Prob) ) )
      END IF
! ----------------------------------------------------------------------
! For the new cloud and precipitation scheme only use the liquid content
! 1. Calculate total of water vapour and liquid water contents, P and
! limit aerosol
! ----------------------------------------------------------------------
      DO Point=1,P_FIELD
        QT(Point) = Q(Point)+QCL(Point)
      END DO

!     Calculate pressure
      DO Point=1,P_FIELD
        P(Point)=P_layer(Point)
      ENDDO

! DEPENDS ON: qsat_wat
      CALL QSAT_WAT (QS,T,P,P_FIELD)

      DO Point = 1 , P_FIELD

!-------------------------------------------------------------------
!* 2. Calculate the ratio of the aerosol mass mixing ratio to the
!*    standard mass mixing ratio, m_over_m0, and the aerosol number
!*    density, N, the dry radius, radius_dry:
!*                      p
!*                  (m )
!*           r = r0 (--)
!*            d     (m0)
!*
!*
!*    And the activation radius:
!*
!*                             1/2
!*                  (       3 )
!*                  ( 3 B0 r  )
!*           r    = ( ------d-)
!*            act   (   A0    )
!*
!*    and A (A0 divided by the dry radius).
! N.B. AEROSOL is in ug/kg, m in kg/kg
! If not available, use 10 ug/kg
!-------------------------------------------------------------------

        if (L_MURK) then
! Ensure that assumed aerosol conc. is at least Aero0:
          m_over_m0 = max(Aerosol(Point)/m0*1.0E-9,                     &
     &                             Aero0/m0*1.0E-9, 0.0001)
        else
          m_over_m0 = max(10.0/m0*1.0E-9, 0.0001)
        endif

        N = N0 * m_over_m0**(1.0-3*power)

        radius_dry = radius0 * (m_over_m0)**power
        A = A0 / radius_dry

        radius_act = SQRT( (3 * B0 * radius_dry**3) / A0 )
        radius_star_act =  radius_act/radius_dry

!-------------------------------------------------------------
!* 3. Calculate the width of the total water
!*    distribution and a modified value of total water, based
!*    on a probability.
!-------------------------------------------------------------

        bs = (1.0-RHcrit) * qs(Point)

        qt_mod = MAX( qt_limit, qt(Point)+ qt_mod_factor* bs)


!====================================================================
!* 4.  Use Newton-Raphson to iteratively improve on a first-guess
!*     droplet radius, using the droplet growth equation and the
!*     geometric relation between liquid water and droplet radius.
!====================================================================
!* 4.1 Calculate a first guess relative humidity, qt/qs, but limit it
!*     to be in the range 0.001 -> 0.999.
!*     From this calculate a first-guess normalised radius using a
!*     simplified version of the droplet growth equation:
!*
!*                              1/3
!*                (       B0   )
!*           r  = ( 1 - ------ )
!*            *   (     ln(RH) )
!*
!----------------------------------------------------------------------

        RH_lim = MIN( MAX( qt_mod/qs(Point), RHmin ) , RHmax )
        radius_star2 = (1.0-B0/LOG(RH_lim))**OneThird

!----------------------------------------------------------------------
!* 4.2 Initialise the iteration counter, the normalised radius
!*     difference, and the updated normalised radius value.
!----------------------------------------------------------------------

        Iteration = 0
        radius_star_diff = 1.0
        radius_star1 = radius_star2

        Do While ( Iteration  <   Niterations .AND.                     &
     &               radius_star_diff  >   Delta_radius_star )

!----------------------------------------------------------------------
!* 4.3 Update the iteration counter and the normalised radius value.
!----------------------------------------------------------------------

          Iteration = Iteration + 1
          radius_star1 = Weight * radius_star2                          &
     &                 + ( 1.0 - Weight ) * radius_star1

!----------------------------------------------------------------------
!* 4.4 Calculate the relative humidity term:
!*
!*                      ( A        B0   )
!*          RHterm = exp( --  -  ------ )
!*                      ( r       3     )
!*                      (  *     r  - 1 )
!*                      (         *     )
!*
!*      and its derivative with respect to the normalised radius:
!*
!*                    (                 2    )
!*                    (   A       3 B0 r     )
!*          RHderiv = ( - --  +  -------*- 2 ) * RHterm
!*                    (    2     (  3     )  )
!*                    (   r      ( r  - 1 )  )
!*                    (    *     (  *     )  )
!*
!----------------------------------------------------------------------

          If ( radius_star1  <   radius_star_act ) then
            RHterm  = EXP( A/radius_star1                               &
     &                     - B0/(radius_star1**3-1.0) )* qs(Point)
            RHderiv = - RHterm * ( -A/(radius_star1**2)                 &
     &                + (3.0*B0*radius_star1**2)                        &
     &                /(radius_star1**3-1.0)**2 )
          Else
            RHterm  = EXP( A/radius_star_act                            &
     &                     - B0/(radius_star_act**3-1.0) ) * qs(Point)
            RHderiv = 0.0
          Endif


!----------------------------------------------------------------------
!* 4.5 Calculate the liquid water mixing ratio term:
!*
!*
!*                   4             3 (  3     )
!*          qLterm = - Pi rho_w N r  ( r  - 1 )
!*                   3             d (  *     )
!*
!*      and its derivative with respect to the normalised radius:
!*
!*                                  3  2
!*          qLderiv = 4 Pi rho_w N r  r
!*                                  d  *
!*
!----------------------------------------------------------------------

          qLterm  = N * FourThirds * Pi * RHO_WATER * radius_dry**3     &
     &                * ( radius_star1**3 - 1.0 )
          qLderiv  = - N * 4.0 * Pi * RHO_WATER                         &
     &                * radius_dry**3 * radius_star1**2

!----------------------------------------------------------------------
!* 4.6 Calculate the function, Fn, and its derivative, Deriv, and
!*     an improved estimate of the normalised radius,
!*     using Newton Raphson:
!*
!*          Fn = qt - RHterm - qLterm
!*
!*          Deriv = RHderiv + qLderiv
!*
!*                          Fn
!*          r      = r  -  -----
!*           * new    *    Deriv
!*
!*     The new estimate of the normalised radius is limited lie between
!*     prescribed maximum and minimum values and within a factor of the
!*     previous value to ensure that the soultion does not diverge.
!----------------------------------------------------------------------

          Fn    = qt_mod - RHterm - qLterm
          Deriv = RHderiv + qLderiv

          radius_star2 = radius_star1 - Fn/Deriv

          IF ( radius_star2  <   radius_star_min )                      &
     &        radius_star2 = radius_star_min
          IF ( radius_star2  >   radius_star_max )                      &
     &        radius_star2 = radius_star_max
          IF ( radius_star2  >   radius_star_factor * radius_star1 )    &
     &        radius_star2 = radius_star_factor * radius_star1
          IF ( radius_star2  <   radius_star1 / radius_star_factor )    &
     &        radius_star2 = radius_star1 / radius_star_factor

!---------------------------------------------------------------------
!* 4.7 Calculate difference between the old and the new values of the
!*     normalised radius.
!---------------------------------------------------------------------

          radius_star_diff = ABS( radius_star1 - radius_star2 )

        END DO

!---------------------------------------------------------------------
!* 5.  Calculate the radius from the final normalised radius.
!---------------------------------------------------------------------

        radius = radius_star2 * radius_dry

!---------------------------------------------------------------------
!* 6. Calculate the visibility, Vis, using the equation:
!*
!*                 ln(liminal contrast)
!*           Vis = -------------2------
!*                     Beta0 N r
!*
!*    (An extra term RecipVisAir is included in the recipical of
!*     visibility to limit visibilities to 100km in clean air).
!---------------------------------------------------------------------

        RecipVis = (N * radius**2) / VisFactor + RecipVisAir
        Visibility(Point) = 1/RecipVis

      END DO

      RETURN
      END SUBROUTINE VISBTY
!
!LL  SUBROUTINE VISTOQT -----------------------------------------------
!LL
!LL     PURPOSE:
!LL Invert relationship between aerosol, visibility and water
!LL content
!LL This is needed for fog probability calculation.
!             Since 4.5, adopted NIMROD based code:
!LL
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.4  07/06/94  First written. Programmer Pete Clark.
!LL  4.0 05/09/95  Diagnosed equivalent RH constrained to be >1%.
!LL                Lower limit to aerosol introduced.  Pete Clark.
!    4.5  30/04/98 NIMROD code adopted.
!                   Pete Clark responsible for UM implementation of
!                   Bruce Wright's NIMROD code.
!    5.1  14/02/00 This routine not activated at 5.1, but change to
!                  aerosol lower bound to be consistent with
!                  visbty routine change. R Rawlins
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3,
!LL                        Version 7, dated 11/3/93.
!LL
!LL  Logical components covered :
!LL
!LL  Project task:
!LL
!LL  External documentation
!LL    Forecasting Research Scientific Paper NO.4
!LL    Diagnosis of visibility in the UK Met Office Mesoscale Model
!LL    and the use of a visibility analysis to constrain initial
!LL    conditions.  SP Ballard, BJ Wright, BW Golding    1992
!      Wright, B. J., 1997: Improvements to the Nimrod Visibility
!         Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!      Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!         for Nimrod. Met. Office FR Tech Rep., No. 222.
!LL
!LLEND----------------------------------------------------------------
!
!*L  Arguments:-------------------------------------------------------
!
!LL  SUBROUTINE FOG_FR------------------------------------------------
!LL
!LL  Purpose: Calculates fog fraction, using the large scale cloud
!LL           scheme. The fog fraction is similar to the cloud fraction
!LL           except it records the fraction of a grid box with RH
!LL           greater than that required for the critical visibility
!LL           (e.g 1 km).
!LL           Since 4.5, adopted NIMROD based code:
!LL           Calculates the fraction of a gridsquare with visibility
!LL           less than threshold, Vis_thresh, given the total water
!LL           mixing ratio, qt, temperature, T, pressure p, and the
!LL           aerosol mass mixing ratio, m, assuming a triangular
!LL           distribution of states about the median, characterised by
!LL           a critical relative humdity value, RHcrit.
!LL           NB:  Throughout, levels are counted from the bottom up,
!LL           i.e. the lowest level under consideration is level 1, the
!LL           next lowest level 2, and so on.
!LL
!LL           Suitable for single-column use.
!LL
!LL   Model         Modification history:
!LL  version  Date
!LL
!LL     3.2 04/05/93 Created by Pete Clark
!LL     3.4 04/08/95 LS_CLD replaced by GLUE_CLD. Andrew Bushell.
!LL   4.2    Oct. 96  T3E migration: *DEF CRAY removed (dynamic
!LL                    allocation now unconditional)
!LL                                   S.J.Swarbrick
!       4.4 01/07/97 Calculation is now based on liquid water
!                    and not on ice. Calculates liquid fog fraction.
!                    This scheme is severely tied to the ideas behind
!                    the 1A cloud scheme. Damian Wilson.
!LL     4.5 30/04/98 NIMROD code adopted. Calculation is still based
!LL                  on liquid water and not ice, but arguments changed
!LL                  to pass T, q, qcl and qcf separately. This is to
!LL                  make code independent of cloud scheme and capable
!LL                  of future development to include ice properly.
!LL                  Pete Clark responsible for UM implementation of
!LL                  Bruce Wright's NIMROD code.
!       5.2 11/10/00 Removed input of num points in global field.
!                    Removed unnecessary error checking.
!                    pstar, ak, bk, replaced by p_layer. D Matthews
!LL
!LL Programming standard:  Unified Model Documentation Paper No 3,
!LL                        Version 5, dated 08/12/92.
!LL
!LL Documentation:
!LL    Wright, B. J., 1997: Improvements to the Nimrod Visibility
!LL       Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!LL    Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!LL       for Nimrod. Met. Office FR Tech Rep., No. 222.
!LL
!LLEND----------------------------------------------------------------
!
!*L
!*LArguments:---------------------------------------------------------
