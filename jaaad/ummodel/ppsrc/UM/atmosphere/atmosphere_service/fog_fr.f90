
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
      SUBROUTINE FOG_FR(                                                &
     & p_layer,RHCRIT,LEVELS,PFIELD,                                    &
     & T,AEROSOL,L_MURK,Q,QCL,QCF,VIS,FF,NVIS)

      IMPLICIT NONE
      INTEGER                                                           &
     & LEVELS                                                           &
                           ! IN No. of levels being processed.
     &,PFIELD                                                           &
                           ! IN No. of points in field (at one level).
     &,NVIS                ! IN No. of visibility thresholds
      REAL                                                              &
     & p_layer(PFIELD,LEVELS)                                           &
                              ! IN pressure (Pa) at processed levels.
     &,RHCRIT(LEVELS)                                                   &
                           ! IN Critical relative humidity.  See the
!                          !    the paragraph incorporating eqs P292.11
!                          !    to P292.14; the values need to be tuned
!                          !    for the given set of levels.
     &,Q(PFIELD,LEVELS)                                                 &
                           ! IN Specific Humidity
!                          !    (kg per kg air).
     &,QCL(PFIELD,LEVELS)                                               &
                           ! Cloud liquid water content at
!                          !     processed levels (kg per kg air).
     &,QCF(PFIELD,LEVELS)                                               &
                           ! Cloud ice content at processed levels
!                          !    (kg per kg air).
     &,T(PFIELD,LEVELS)                                                 &
                           ! IN Temperature (K).
     &,AEROSOL(PFIELD,LEVELS)                                           &
                              ! IN Aerosol mixing ratio(ug/kg)
     &,VIS(PFIELD,LEVELS,NVIS)  ! Visibility thresholds
      LOGICAL                                                           &
     &   L_MURK               ! IN : Aerosol present

      REAL                                                              &
     & FF(PFIELD,LEVELS,NVIS)   ! OUT Vis prob at processed levels
!                          !     (decimal fraction).
!
!*--------------------------------------------------------------------
!*L  Workspace usage----------------------------------------------------
      REAL                                                              &
                           ! "Automatic" arrays on Cray.
     & P(PFIELD)                                                        &
     &,QT(PFIELD)                                                       &
                           ! total of cloud water and vapour
     &,QS(PFIELD)                                                       &
                           ! Saturated spec humidity for temp T
     &,qt_thresh(PFIELD)                                                &
                           ! modified qt
     &,bs
!*L  External subroutine called ----------------------------------------
      EXTERNAL QSAT_WAT,VISTOQT
!* Local, including SAVE'd, storage------------------------------------
!
      INTEGER K,I,J     ! Loop counters: K - vertical level index.
!                       !                I - horizontal field index.
                        !                J - Vis threshold index.
!*--------------------------------------------------------------------
!*  Local and other physical constants----------------------------------
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
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
!
!-----------------------------------------------------------------------
!L Subroutine structure :
!L Loop round levels to be processed.
!-----------------------------------------------------------------------
!
      DO K=1,LEVELS
!
!-----------------------------------------------------------------------
!L 1. Calculate Pressure and initialise temporary arrays
!-----------------------------------------------------------------------
!
        DO I=1,PFIELD
          P(I)=p_layer(I,K)
          QT(I)=Q(I,K)+QCL(I,K)
        ENDDO ! Loop over points

!-----------------------------------------------------------------------
!* 2.  Calculate total water threshold corresponding to visibility
!      Since Qs is needed more than once, pre-calculate and pass it
!-----------------------------------------------------------------------

! DEPENDS ON: qsat_wat
        CALL QSAT_WAT (Qs,T(1,K),P,PFIELD)

        DO J=1,NVIS

! DEPENDS ON: vistoqt
          Call VISTOQT( VIS(1,K,J), Qs, AEROSOL(1,K), L_MURK,           &
     &                PFIELD, qt_thresh )


!-----------------------------------------------------------------------
!* 3.  Calculate the width of the distribution in total water space, bs:
!*
!*           bs = ( 1 - RHcrit ) * qs(T)
!*
!-----------------------------------------------------------------------

          Do I = 1 , PFIELD

            bs = (1.0-RHcrit(K)) * qs(I)

!=======================================================================
!* 4.  Calculate the fraction of states in a triangular
!*     distribution which exceed the total water threshold.
!=======================================================================

!-----------------------------------------------------------------------
!* 4.1 If total water threshold value is less than the total water value
!*     minus the width of the distribution, then all of the states have
!*     a total water value exceeding the threshold, so set the
!*     visibility fraction to 1.0
!-----------------------------------------------------------------------

            if ( qt_thresh(I)  <=  qt(I)-bs ) then

              FF(I,K,J) = 1.0

!-----------------------------------------------------------------------
!* 4.2 If total water threshold value is greater than the total water
!*     value minus the width of the distribution, but less than the
!*     total water value then the visibility fraction, VF, is given by:
!*
!*                                                    2
!*                             ( qt       - qt + bs  )
!*            VF = 1.0 - 0.5 * (    thresh           )
!*                             ( ------------------- )
!*                             (          bs         )
!*
!-----------------------------------------------------------------------

             Else if ( qt_thresh(I)  >   qt(I)-bs .AND.                 &
     &                 qt_thresh(I)  <=  qt(I) ) then

               FF(I,K,J) = 1.0 - 0.5 *                                  &
     &              (( qt_thresh(I) - qt(I) + bs )/ bs)**2

!-----------------------------------------------------------------------
!* 4.3 If total water threshold value is greater than the total water
!*     value, but less than the total water value plus the width of the
!*     distribution, then the visibility fraction, VF, is given by:
!*
!*                                              2
!*                       ( qt + bs - qt        )
!*            VF = 0.5 * (             thresh  )
!*                       ( ------------------- )
!*                       (          bs         )
!*
!-----------------------------------------------------------------------

             Else if ( qt_thresh(I)  >   qt(I) .AND.                    &
     &                 qt_thresh(I)  <=  qt(I)+bs    ) then

                FF(I,K,J)= 0.5 * (( qt(I) + bs - qt_thresh(I))/bs)**2

!-----------------------------------------------------------------------
!* 4.4 If total water threshold value is greater than the total water
!*     value plus the width of the distribution, then non of the states
!*     have a total water value exceeding the threshold, so set the
!*     visibility fraction to 0.0
!-----------------------------------------------------------------------

             Else

               FF(I,K,J) = 0.0

            End if

          End Do ! Loop over PFIELD I

        End Do ! Loop over VIS J

      ENDDO ! Loop over levels
!
      RETURN
      END SUBROUTINE FOG_FR
