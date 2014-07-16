#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate specific humidities from Gill's formula.
!
! Method:
!       Straightforward.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE QSAT_GILL(SAT_SPEC_HUM, T, P                           &
     &  , N_PROFILE, N_LAYER                                            &
     &  , ND_PROFILE, ND_LAYER                                          &
     &  )
!
!
      IMPLICIT NONE
!
! Include Header Files
#include "c_kinds.h"
!
!     This routine computes the saturated specific humidity
!     at temperature T and pressure p, using the formulae given in
!     Appendix four of Adrian Gill's book.
!
!     Note that the formulae work with pressures in hectopascals and
!     temperatures in degrees celsius. the conversions are made
!     inside this routine and should have no impact on the rest
!     of the code.
!
!     This routine was perpetrated by D L Roberts (12/8/93).
!
!     Modified to cope with very low pressures by DLR (27/10/93).
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles.
     &  , N_LAYER                                                       &
!           Number of layers.
     &  , ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles.
     &  , ND_LAYER
!           Size allocated for atmospheric layers.
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure in pascals
     &  , T(ND_PROFILE, ND_LAYER)
!           Temperature in kelvin
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    SAT_SPEC_HUM(ND_PROFILE, ND_LAYER)
!           Saturated specific humidity at T and p.
!
!
!     Local variables.
!
      INTEGER I,L
!
      REAL  (Real64) :: PRESS
!       PRESSURE IN HECTOPASCALS.
      REAL  (Real64) :: TEMP
!       TEMPERATURE IN CELSIUS.
      REAL  (Real64) :: A
!       A TEMPORARY HOLDING VARIABLE.
      REAL  (Real64) :: EW
!       SATURATION VAPOUR PRESSURE OF PURE WATER VAPOUR.
      REAL  (Real64) :: EWDASH
!       SAT VAP PRESSURE OF WATER VAPOUR IN AIR.
      REAL  (Real64) :: FW
!       THE RATIO BETWEEN EWDASH AND EW.
      REAL  (Real64) :: ZERO_DEGC
!       KELVIN EQUIVALENT OF ZERO CELSIUS.
!
!     The value assigned is that used in v3.1 of the unified model.
!
      PARAMETER( ZERO_DEGC=273.15E+00_Real64 )
!
      REAL  (Real64) :: EPSILON ! THE RATIO OF THE
!                              ! MOLECULAR MASS OF WATER
!                            to that of dry air.
!     The value assigned is that used in v3.1 of the unified model.
!
      PARAMETER( EPSILON=0.62198E+00_Real64 )
!
      REAL  (Real64) :: ETA ! ONE MINUS EPSILON.
      PARAMETER( ETA = 1.0E+00_Real64 - EPSILON )
!
!
!     Loop over all points.
!     These loops are not indented, in order to make the
!     equations easier to read by keeping them to one line.
!
      DO I = 1,N_LAYER
      DO L = 1,N_PROFILE
!
!     Convert to local units for temperature and pressure.
!
      TEMP = T(L,I) - ZERO_DEGC
      PRESS = P(L,I)*0.01E+00_Real64
!
!     Equation (A4.7) of Gill's book.
!
      FW = 1.0E+00_Real64 + 1.0E-06_Real64*PRESS*( 4.5E+00_Real64       &
     &  + 6.0E-04_Real64*TEMP*TEMP )
!
!     Equation (A4.5) of Gill.
!
      A = ( 7.859E-01_Real64 + 3.477E-02_Real64*TEMP )                  &
     &   /( 1.0E+00_Real64 + 4.12E-03_Real64*TEMP )
      EW = 1.0E+01_Real64**A
!
!     Equation (A4.6) of Gill.
!
      EWDASH = FW*EW
!
!     The next equation is a rearrangement of Gill's (A4.3),
!     with w subscripts added because we require saturation.
!
!     Note that at very low pressures a fix has to be applied,
!     to avoid a singularity.
!
      IF (PRESS  >   EWDASH) THEN
        SAT_SPEC_HUM(L,I) = (EPSILON*EWDASH)/(PRESS-EWDASH*ETA)
      ELSE
        SAT_SPEC_HUM(L,I) = 1.0E+00_Real64
      ENDIF
!
!     End of the double loop.
!
      ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE QSAT_GILL
#endif
