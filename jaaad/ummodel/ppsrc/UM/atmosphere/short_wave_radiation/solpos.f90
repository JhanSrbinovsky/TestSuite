
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
      SUBROUTINE SOLPOS (DAY, YEAR, LCAL360, L_SEC_VAR, L_EQT           &
     &                 , EQT, SINDEC, SCS )
!
      IMPLICIT NONE
!
      LOGICAL, INTENT(IN) :: LCAL360    !True if 360 day calendar in use
      LOGICAL, INTENT(IN) :: L_EQT      !True, include equation of time
      LOGICAL, INTENT(IN) :: L_SEC_VAR  !True, include secular variation
!                                       !      of the Earth's orbit
      INTEGER, INTENT(IN) :: DAY        ! Day-number in the year
      INTEGER, INTENT(IN) :: YEAR       ! Calendar year
!
      REAL, INTENT(OUT)   :: SINDEC     ! Sin(solar declination)
      REAL, INTENT(OUT)   :: SCS        ! Solar constant scaling factor
      REAL, INTENT(OUT)   :: EQT        ! The equation of time,
!                                       !  specified as an hour angle
!                                       !  in radians.

!L This routine has no dynamically allocated work areas and no
!L  significant structure.  It calls the intrinsic functions FLOAT, SIN
!   & COS, and the user subroutine ORBPRM
!  The user deck <astron/astron.h> is included.
!
!
!     Mathematical constants:
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
      REAL, PARAMETER :: TWOPI = 2. * PI
!
!     Parameters of the Earth's orbit:
      REAL :: E                        ! Eccentricity of the orbit
      REAL :: GAMMA                    ! Supplement of the longitude
!                                      !  of the perihelion
      REAL :: TAU0                     ! Time of the perihelion
!                                      !  passage in days
      REAL :: OBLQ                     ! Obliquity of the orbit
      REAL :: DINY                     ! Length of the calendar year
!                                      !  (in whole days)
!
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! ----------------------- Header file YEARLEN  -------------------------
! Description: Parameter of the length of the tropical year
!
! Current Code Owner: J. M. Edwards
!
! History:
! Version  Date      Comment.
!  5.2     25/10/00  Original Code.   E. Ostrom
!
!----------------------------------------------------------------------
!     ! Number of days in the tropical year
!     ! The tropical year is defined as the mean interval between two
!     ! successive passages of the sun through the vernal equinox.
!     ! This value is 365.2424 and will be used here.
!     ! The value 365.2422 most often used is the mean length of
!     ! the year starting at different points on the ellipse.
      Real, Parameter    :: TropYearLength = 365.2424
!
!     Derived orbital constants:
      REAL :: E1
      REAL :: E2
      REAL :: E3
      REAL :: E4
      REAL :: Y
!
      REAL :: M       ! Mean anomaly: positional angle of a "mean" Earth
!                     !  rotating around the sun with a constant angular
!                     !  speed equal to 2pi/T and counted counterclock-
!                     !  wise from the perihelion
      REAL :: V       ! True anomaly: positional angle of Earth in its
!                     !  orbit, counted countercloskwise from perihelion
!
!     Determine the orbital parameters for the calendar selected.
! DEPENDS ON: orbprm
      CALL ORBPRM(L_SEC_VAR, YEAR, LCAL360                              &
     &          , E, GAMMA, OBLQ, TAU0, DINY )
!
!     Calculate the mean anomaly at 12Z on the current day.
!     The 0.5 accounts for the time in days since mid-night.
!     The references are to Smart 1944 (and UMDP23)
!     Eq 67 p. 113 and n=2pi/orbital period     (Eq 3.1.1)

      IF (LCAL360) THEN
        M = (TWOPI / DINY)           * (FLOAT(DAY) - TAU0 - .5)
      ELSE
        M = (TWOPI / TropYearLength) * (FLOAT(DAY) - TAU0 - .5)
      END IF

!       Calculate the coefficients in the equation of the centre and
!        thence derive the true anomaly.
      E1 = E * ( 2. - .25 * E*E )
      E2 = 1.25 * E*E
      E3 = E*E*E * 13. / 12.

!       True anomaly, equation 87 in Smart on p. 120 (UMDP23 Eq 3.1.2)
      V  = M + E1*SIN(M) + E2*SIN(2.*M) + E3*SIN(3.*M)

!       Solar constant scaling factor (UMDP23 Eq 3.1.4)
      E4  = ( (1. + E*E*.5) / (1. - E*E) )**2
      SCS = E4 * ( 1. + E * COS(V) ) **2

!       sin(solar declination) (UMDP23 Eq 3.1.5)
!       The solar declination is related to
!        the true longitude of the earth (lambda) by:
!        sindec = sin(obliquity) * sin(lambda)
!       Lambda is counted counterclockwise from the vernal equinox
!        and is related to v (the true anomaly) through
!        lambda = v + (longitude of perihelion)

      SINDEC = SIN(OBLQ) * SIN (V - GAMMA)
!
!
!     Calculate the equation of time as given by equation (29)
!     on page 149 of Smart (1944). (Recall the factor of 5/4 in the
!     definition of E2).
      IF (L_EQT) THEN
        Y   = ( TAN ( 0.5*OBLQ ) )**2
        EQT = Y * SIN(2.0 * ( M - GAMMA ))                              &
     &        - 2.0*E * SIN(M)                                          &
     &        + 4.0*E*Y * SIN(M) * COS(2.0*( M - GAMMA ))               &
     &        - 0.5*Y*Y * SIN(4.0*( M - GAMMA ))                        &
     &        - E2 * SIN(2.0*M)
      ELSE
        EQT=0.0E+00
      END IF

      RETURN
      END SUBROUTINE SOLPOS
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
