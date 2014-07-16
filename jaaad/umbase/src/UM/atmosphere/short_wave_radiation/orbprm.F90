#if defined(A01_3A)||defined(A01_3C)||defined(A01_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine SOLPOS   ----------------------------------------------
!LL
!LL Purpose :
!LL  Calculations of the earth's orbit described in the first page of
!LL  the "Calculation of incoming insolation" section of UMDP 23, i.e.
!LL  from the day of the year (and, in forecast mode, whether it is a
!LL  leap year) and the orbital "constants" (which vary over
!LL  "Milankovitch" timescales) it calculates the sin of the solar
!LL  declination and the inverse-square scaling factor for the solar
!LL  "constant".  It is thus intrinsically scalar.  The FORTRAN code
!LL  present depends on whether *DEF CAL360 is set during UPDATE: this
!LL  replaces the Julian calendar with the climate-mode 360-day calendar
!LL
!LL   Author:    William Ingram
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.4    20/06/94 DEF CAL360 replaced by LOGICAL LCAL360;
!LL                   PARAMETER statements duplicated for 360 and
!LL                   365 day calendar.
!LL                                                S.J.Swarbrick
!LL   4.4    27/02/97 Testing for leap years modified to deal with
!LL                   no leap every 100y except for every 400y
!LL                   Author: M.Gallani
!LL
!LL   5.1    14/06/99 Modified to include the equation of time and
!LL                   more modern values of astronomical constants
!LL                   as recommended for AMIP-II. The equation of
!LL                   time (the difference netween true and mean
!LL                   solar time) is calculated once for the day.
!LL                                                J. M. Edwards
!
!     5.2    15/11/00 Modified to include the secular variation of
!                     the orbital parameters, if that option is
!                     chosen. Use the logical switch L_SEC_VAR.
!                                                  E. Ostrom
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
!   Programming standard:
!      Written in FORTRAN90 to comply with Version 7.2 of the UMDP3
!       dated the 5/2/98.
!LL Logical components covered : P233
!LL
!LL Project task :
!LL
!LL External documentation: P23
!LL
!LLEND -----------------------------------------------------------------
!*L
!
!+ Subroutine to calculate the parameters of the Earth's orbit.
!
!  Purpose:
!  This routine returns the parameters of the Earth's orbit
!   (the eccentricity, obliquity, supplement of the longitude of
!    perihelion) the time of the parihelion passage in days, and
!   Length of the calendar year (in whole days).
!
!  Method:
!  For long runs there may be an interest in running with secular
!   variations in the astronomy. The orbital constants have
!   been derived from A. L. Berger 1978, J. Atm. Sci, Volume 35
!   2362-2367. A copy of which can be found in the Met Office library.
!  For short current runs, or long control runs it is preferrable
!   not to allow the astronomy to vary, so fixed values are used.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       5.2             15/11/00                Original Code
!                                               E. Ostrom
!
! Description of Code:
!   FORTRAN90 complying with UMDP3 Version 7.2 from 5/2/98
!
!- ---------------------------------------------------------------------
      SUBROUTINE ORBPRM(L_SEC_VAR, YEAR, LCAL360                        &
     &                , E, GAMMA, OBLQ, TAU0, DINY )

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: YEAR       ! Calendar year
      LOGICAL, INTENT (IN) :: L_SEC_VAR  ! Include secular variations
!                                        !  of the orbit
      LOGICAL, INTENT (IN) :: LCAL360    ! Use a calendar of 360 days

!     Parameters of the Earth's orbit:
!
      REAL, INTENT(OUT) :: E             ! Eccentricity of the orbit
      REAL, INTENT(OUT) :: GAMMA         ! Supplement of the longitude
!                                        !  of the perihelion
      REAL, INTENT(OUT) :: OBLQ          ! Obliquity of the orbit
      REAL, INTENT(OUT) :: TAU0          ! Time of the perihelion
!                                        !  passage in days
      REAL, INTENT(OUT) :: DINY          ! Length of the calendar year
!                                        !  (in whole days)

!     Local Variables for use within ORBPRM
!
      REAL :: YEAR_OFFSET                ! Offset of the year from the
!                                        !  reference year when default
!                                        !  values apply
      REAL :: ECN_SN                     ! Eccentricity multiplied by
!                                        !  the sine of the longitude
!                                        !  of the perihelion
      REAL :: ECN_CN                     ! Eccentricity multiplied by
!                                        !  the cosine of the longitude
!                                        !  of the perihelion
      REAL :: LPH_FIXED_VE               ! Longitude of the perihelion
!                                        !  relative to a fixed vernal
!                                        !  equinox
      REAL :: GN_PRCS                    ! General precession
      REAL :: DATE_VE                    ! Date of the vernal equinox
!                                        !  in days into the year
      REAL :: NO_LEAP_DAYS               ! The number of leap days,
!                                        !  used to calculate DATE_VE
      REAL :: MEAN_ANOM_VE               ! Mean anomaly at the vernal
!                                        !  equinox

!     Synthetic constants
!
      REAL :: BETA
      REAL :: EE1
      REAL :: EE2
      REAL :: EE3

      INTEGER :: I                       ! Loop variable

!     Mathematical constants:
#include "c_pi.h"
      REAL, PARAMETER :: TWOPI = 2. * PI

!     Astronomical Parameters:
#include "astron.h"
#include "yearlen.h"

!     The length of the calendar year may be set for a 360-day calendar
!      (as is often used in climate runs),
!      or for a real Gregorian calendar which has 365 days in
!      non-leap years and 366 in leap years.

      IF (LCAL360) THEN

        DINY=360.0

      ELSE
!      Is this a leap year?
        IF (mod(year,4)    ==  0 .AND.                                  &
     &     (mod(year,400)  ==  0 .OR. mod(year,100)  /=  0)) then

          DINY = 366.0

!      Is this a normal year?
        ELSE

          DINY = 365.0

        END IF
      END IF

!     The orbital elements are normally set to default values, but
!     secular variations may be required in some longer climate runs.

      IF (L_SEC_VAR) THEN

        YEAR_OFFSET = REAL( YEAR - YEAR_REF )

!       Obliquity: (Equation 1 from Berger 1978)

        OBLQ = OBLQ_CNST
        DO I=1, N_TERM_OBQ
          OBLQ = OBLQ+A(I)*COS(F(I)*YEAR_OFFSET+D(I))
        END DO

!       Eccentricity: this is better computed from its components
!       than directly from the series.(Equation (4) of Berger 1978).

        ECN_SN = M(1) * SIN (G(1) * YEAR_OFFSET + B(1))
        ECN_CN = M(1) * COS (G(1) * YEAR_OFFSET + B(1))

        DO I=2, N_TERM_ECN_LPH
          ECN_SN = ECN_SN + M(I) * SIN (G(I) * YEAR_OFFSET + B(I))
          ECN_CN = ECN_CN + M(I) * COS (G(I) * YEAR_OFFSET + B(I))
        END DO
        E = SQRT(ECN_SN*ECN_SN+ECN_CN*ECN_CN)

!       We now obtain the longitude of the perihelion relative to the
!       fixed equinox.

        LPH_FIXED_VE = ATAN2 (ECN_SN,ECN_CN)

!       The longitude of perihelion and
!        the supplement of the longitude of the perihelion relative to
!        the actual vernal equinox requires the general precession.

!      General Precession.
        GN_PRCS = LIN_RATE_GN_PRCS * YEAR_OFFSET + GN_PRCS_CNST
        DO I=1, N_TERM_GN_PRCS
          GN_PRCS = GN_PRCS + C(I) * SIN (H(I) * YEAR_OFFSET + R(I))
        END DO

!      Supplement of the longitude of the perihelion
        GAMMA = PI - LPH_FIXED_VE - GN_PRCS

!       Time of perihelion: The time at which an object is at perihelion
!        (its closest distance to the sun).
!       The time of perihelion is inferred from the date of
!        the vernal equinox using the Gregorian calendar.
!
!      Calculate the date of the vernal equinox.
!       First we need to:
!        Calculate the no of leap days between year & year_ref_ve.
!        This needs to be corrected when using the Gregorian calendar.
!         by adding (DINY-366.0) when the year_ref_ve is a leap year or
!         by adding (DINY-365.0) when the year_ref_ve is a normal year.
!        This correction is done when the DATE_VE is calculated below!
!
!        In the calculation of NO_LEAP_DAYS below, the divisions of type
!         'YEAR'/x (where x is 4, 100 or 400) are integer computations.
!         These integers are then subtracted and the resulting integer
!         is then converted to a real.

        NO_LEAP_DAYS = ( TropYearLength - 365.0)                        &
     &    * REAL( YEAR     - YEAR_REF_VE     )                          &
     &    - REAL( YEAR/4   - YEAR_REF_VE/4   )                          &
     &    + REAL( YEAR/100 - YEAR_REF_VE/100 )                          &
     &    - REAL( YEAR/400 - YEAR_REF_VE/400 )

!      Now we can calculate the date of the vernal equinox!
!      Because the date of the vernal equinox is varying with the year,
!      we have to keep track of its position in the sky.
!      In order to accomodate a time varying vernal equinox when using
!      a 360-day year, we still have to calculate the difference in
!      the vernal equinox depending on leap years, normal years and the
!      difference between the length of the tropical year and the
!      "normal" year and then we adjust this by multiplying the
!      DATE_VE by 360/(length of tropical year).
!
!      Is a 360 day calendar being used?

        IF (LCAL360) THEN
          DATE_VE = DATE_VE_DFLT + NO_LEAP_DAYS
          DATE_VE = DATE_VE * DINY / TropYearLength


!      Is a 365 day calendar being used?
        ELSE

!        Is the epoch reference year a leap year?

          IF (mod(YEAR_REF_VE,4)    ==  0 .AND.                         &
     &       (mod(YEAR_REF_VE,400)  ==  0 .OR.                          &
     &        mod(YEAR_REF_VE,100)  /=  0)) THEN

            DATE_VE = DATE_VE_DFLT + (NO_LEAP_DAYS + (DINY - 366.0))

!        Is the epoch reference year a normal year?

          ELSE

            DATE_VE = DATE_VE_DFLT + (NO_LEAP_DAYS + (DINY - 365.0))

          END IF
        END IF

        BETA = SQRT(1.0E+00-E*E)
        EE1  = (0.5*E + 0.125*E*E*E)*(1.0 + BETA)
        EE2  = -0.25*E*E* (0.5 + BETA)
        EE3  = 0.125*E*E*E*((1.0/3.0) + BETA)
        MEAN_ANOM_VE = GAMMA - 2.0E+00 * (                              &
     &      EE1 * SIN (GAMMA)                                           &
     &    + EE2 * SIN (2.0 * GAMMA)                                     &
     &    + EE3 * SIN (3.0 * GAMMA)                                     &
     &    )

        TAU0 = DATE_VE - MEAN_ANOM_VE * TropYearLength/(TWOPI)

      ELSE

        E     = E_DFLT
        OBLQ  = OBLQ_DFLT
        GAMMA = PI - LPH_DFLT
        TAU0  = TAU0_DFLT

      END IF

!     If using a 360-day calendar the time of the perihelion is
!     adjusted.
      IF (LCAL360) THEN
        TAU0 = TAU0*(360.0/TropYearLength)+0.71
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE ORBPRM
#endif
