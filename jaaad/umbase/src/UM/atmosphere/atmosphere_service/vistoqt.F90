#if defined(A71_1A)
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
!    6.2  03/02/06  Moved to A71_1A. P.Selwood
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
      SUBROUTINE VISTOQT                                                &
     &           (VISIBILITY                                            &
                                                      !INPUT
     &           ,Qs                                                    &
                                                      !INPUT
     &           ,AEROSOL                                               &
                                                      !INPUT
     &           ,L_MURK                                                &
                                                      !INPUT
     &           ,Npoints                                               &
                                                      !INPUT
     &           ,qt )                                !OUTPUT
      IMPLICIT NONE
!---------------------------------------------------------------------
! Workspace usage:----------------------------------------------------
! None
!*--------------------------------------------------------------------
!*L-------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
      INTEGER                                                           &
     &        Npoints             ! IN NO. points in field.
      REAL                                                              &
     &        VISIBILITY(Npoints)                                       &
                                 ! IN visibility
     &       ,Qs(Npoints)                                               &
                                  !  Saturated humidity mixing ratio
     &       ,AEROSOL(Npoints)    ! IN Aerosol mixing ratio(ug/kg)
      LOGICAL                                                           &
     &        L_MURK                    ! IN : Aerosol present
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
      REAL                                                              &
     &        qt(Npoints)         ! OUT Total water mixing ratio (kg/kg)
!*--------------------------------------------------------------------
!*L-------------------------------------------------------------------
!---------------------------------------------------------------------
!*--------------------------------------------------------------------
! constants for visibility calculation used to be set here but now
! set in COMDECK.
! PI needed to set new constants.
!---------------------------------------------------------------------
#include "c_pi.h"
#include "c_r_cp.h"
#include "c_epslon.h"
#include "c_densty.h"
#include "c_visbty.h"
!
! Local parameter variables
!
      REAL                                                              &
     &  qt_limit       ! Smallest total water mixing ratio value allowed
      PARAMETER (   qt_limit = 0.0001 )
!
! Local workspace variables
!
       INTEGER                                                          &
     &  Point            !  Loop variable for points

      REAL                                                              &
     &  qL                                                              &
                         !  Liquid water mixing ratio (Kg/Kg).
     &, radius_dry                                                      &
                         !  Dry particle radius for aerosol (m)
     &, radius                                                          &
                         !  Radius of fog droplets (m)
     &, radius_star                                                     &
                         !  Water droplet radius divided by dry radius
     &, radius_act                                                      &
                         !  Activation droplet radius
     &, radius_star_act                                                 &
                         !  Activation droplet radius divided by the dry
     &, radius_star_used                                                &
                         !  Water droplet radius divided by the dry
                         !   radius actually used for the relative
                         !   humidity calculation
     &, RH                                                              &
                         !  Relative humidity derived from visibility
     &, A                                                               &
                         !  A0 divided by the dry radius
     &, m_over_m0                                                       &
                         !  Ratio of the aerosol mass mixing ratio and
                         !   the standard aerosol mass mixing ratio
     &, N                !  Number density of aerosol particles (/m3)


!**

      Do Point = 1 , Npoints

!---------------------------------------------------------------------
!* 1.  Calculate the ratio of the aerosol mass mixing ratio to the
!*     standard mass mixing ratio, m_over_m0, and the aerosol number
!*     density, N:
!*
!*                      (1-3p)
!*                  (m )
!*           N = N0 (--)
!*                  (m0)
!*
!*     And the dry radius, radius_dry:
!*                      p
!*                  (m )
!*           r = r0 (--)
!*            d     (m0)
!*
!*     And A (A0 divided by the dry radius).
!*
!*     And the activation radius:
!*
!*                             1/2
!*                  (       3 )
!*                  ( 3 B0 r  )
!*           r    = ( ------d-)
!*            act   (   A0    )
!*
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

!----------------------------------------------------------------------
!* 2.  Calculate a water droplet radius, from the visibility:
!*
!*                                1/2
!*               ( ln( epsilon ) )
!*           r = (---------------)
!*               (  Vis N Beta0  )
!*
!*    (An extra term RecipVisAir is included in the recipical of
!*     visibility to limit visibilities to 100km in clean air).
!----------------------------------------------------------------------

              radius =(VisFactor/N) *                                   &
     &           (1.0/Visibility(Point) - RecipVisAir)
              IF (radius  >   0.0) THEN
                radius = SQRT( radius )
              ELSE
                radius = radius_dry
              ENDIF

!----------------------------------------------------------------------
!* 3.  Provided the diagnosed radius is greater than the dry radius,
!*     calculate the normalised droplet radius, and the saturated
!*     humidity mixing ratio.
!----------------------------------------------------------------------

        If ( radius  >   radius_dry ) then

          radius_star = radius / radius_dry

!----------------------------------------------------------------------
!* 5.  Calculate the corresponding liquid water mixing ratio:
!*
!*
!*               4            (  3     3 )
!*          qL = - Pi rho_w N ( r  - r   )
!*               3            (       d  )
!*
!----------------------------------------------------------------------

          qL = FourThirds * Pi * rho_water * N  *                       &
     &         ( radius**3 - radius_dry**3 )

!----------------------------------------------------------------------
!* 6.  Calculate the relative humidity:
!*
!*                  ( A        B0   )
!*          RH = exp( --  -  ------ )
!*                  ( r       3     )
!*                  (  *     r  - 1 )
!*                  (         *     )
!*
!----------------------------------------------------------------------

          If ( radius_star  <   radius_star_act ) then
            RH = EXP( A/radius_star                                     &
     &                - B0 /( radius_star **3 - 1.0 ) )
          Else
            RH = EXP( A/radius_star_act                                 &
     &                  - B0 /( radius_star_act **3 - 1.0 ) )
          Endif

!----------------------------------------------------------------------
!* 7.  Calculate the total water mixing ratio:  qt = RH * qs(T) + qL
!----------------------------------------------------------------------

          qt(Point) = MAX( RH * qs(Point) + qL, qt_limit )

!----------------------------------------------------------------------
!* 8. If the droplet radius is less than the dry radius, then set the
!*    total water mixing ratio to the minimum value.
!----------------------------------------------------------------------

        Else

          qt(Point) = qt_limit

        End if


      End Do

      RETURN
      END SUBROUTINE VISTOQT
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
#endif
