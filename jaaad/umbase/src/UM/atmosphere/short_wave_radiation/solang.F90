#if defined(A01_3A) || defined(A01_3C) || defined(A01_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculation of solar zenith/hour angles and sunlit fraction.
!
! Purpose :
!  Calculations of the earth's orbit described in the second page of
!  the "Calculation of incoming insolation" section of UMDP 23, i.e.
!  from the sin of the solar  declination, the position of each point
!  and the time limits it calculates how much sunlight, if any, it
!  receives.
!
! Original Author:    William Ingram
!
! Current Code Owner: James Manners
!
!- ---------------------------------------------------------------------
!
      SUBROUTINE SOLANG (SINDEC, T, DT, EQT, SINLAT, LONGIT, K,         &
     &     LIT, COSZ, MEAN_OMEGA)
!
      USE solinc_data, ONLY: sol_bearing
      IMPLICIT NONE
!
      INTEGER                                                           &
             !, INTENT(IN) ::
     &     K                          ! Number of points
      REAL                                                              &
          !, INTENT(IN) ::
     &     SINDEC,                                                      &
                                      ! Sin(solar declination)
     &     T, DT,                                                       &
                                      ! Start time (GMT) & timestep
     &     SINLAT(K),                                                   &
                                      ! sin(latitude) & longitude
     &     LONGIT(K),                                                   &
                                      ! of each point
     &     EQT                        ! The equation of time (the
!                                     ! difference between true and
!                                     ! mean solar time): this is
                                      ! calculated only once each day
                                      ! like the solar declination: it
                                      ! is set to 0 if not to be used.
      REAL                                                              &
          !, INTENT(OUT) ::
     &     LIT(K),                                                      &
                                      ! Sunlit fraction of the timestep
     &     COSZ(K),                                                     &
                                      ! Mean cos(solar zenith angle)
!                                     ! during the sunlit fraction
     &     MEAN_OMEGA(K)
                           ! Mean hour angle over the timestep
                           ! measured in radians west of local noon.
!*
!L This routine has no dynamically allocated work areas.  It calls the
!L intrinsic functions SQRT, ACOS & SIN, but no user functions or
!L subroutines.  The only structure is a loop over all the points to be
!L dealt with, with IF blocks nested inside to cover the various
!L possibilities.
      INTEGER J                       ! Loop counter over points
      REAL TWOPI,                                                       &
                                      ! 2*pi
     &     S2R                        ! Seconds-to-radians converter
      REAL SINSIN,                                                      &
                              ! Products of the sines and of the cosines
     &     COSCOS,                                                      &
                              ! of solar declination and of latitude.
     &     HLD,                                                         &
                              ! Half-length of the day in radians (equal
                              ! to the hour-angle of sunset, and minus
     &     COSHLD,                                                      &
                              ! the hour-angle of sunrise) & its cosine.
     &     HAT,                                                         &
                              ! Local hour angle at the start time.
     &     OMEGAB,                                                      &
                              ! Beginning and end of the timestep and
     &     OMEGAE,                                                      &
                              ! of the period over which cosz is
     &     OMEGA1,                                                      &
                              ! integrated, and sunset - all measured in
     &     OMEGA2,                                                      &
                              ! radians after local sunrise, not from
     &     OMEGAS,                                                      &
                              ! local noon as the true hour angle is.
     &     DIFSIN,                                                      &
                              ! A difference-of-sines intermediate value
     &     DIFTIM,                                                      &
                              ! and the corresponding time period
     &     TRAD, DTRAD                                                  &
!     ! These are the start-time and length of the timestep (T & DT)
!     ! converted to radians after midday GMT, or equivalently, hour
!     ! angle of the mean sun on the Greenwich meridian.
     &     ,X(K)                                                        &
     &     ,LAT(K),DEC(K)  ! Working variables for bearing calculation.
#include "c_pi.h"
      PARAMETER ( TWOPI = 2. * PI, S2R = PI / 43200.)
!
      TRAD = T * S2R - PI
      DTRAD = DT * S2R
!DIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO 100 J = 1, K                          ! Loop over points
       COSHLD = 0.
       HLD = 0.                                ! Logically unnecessary
! statement without which the CRAY compiler will not vectorize this code
       SINSIN = SINDEC * SINLAT(J)
       COSCOS = SQRT( (1.-SINDEC**2) * (1.-SINLAT(J)**2) )
       IF ( SINSIN  <   -COSCOS ) THEN                   ! Perpetual nig
          LIT(J) = 0.
          COSZ(J) = 0.
          MEAN_OMEGA(J) = 0.
        ELSE
!         Since T and DT represent mean time and all solar calculations
!         are done using the variable HAT, it suffices to add the
!         equation of time on to HAT.
          HAT = LONGIT(J) + TRAD + EQT         ! (3.2.2)
          IF ( SINSIN  >   COSCOS ) THEN
             OMEGA1 = HAT                      ! angles for (3.2.3) are
             OMEGA2 = HAT + DTRAD              ! start & end of timestep
           ELSE                                !   At this latitude some
             COSHLD = SINSIN / COSCOS
! points are sunlit, some not.  Different ones need different treatment.
             HLD = ACOS(-COSHLD)               ! (3.2.4)
! The logic seems simplest if one takes all "times" - actually hour
! angles - relative to sunrise (or sunset), but they must be kept in the
! range 0 to 2pi for the tests on their orders to work.
             OMEGAB = HAT + HLD
             IF (OMEGAB <  0.)   OMEGAB = OMEGAB + TWOPI
             IF (OMEGAB >= TWOPI) OMEGAB = OMEGAB - TWOPI
             IF (OMEGAB >= TWOPI) OMEGAB = OMEGAB - TWOPI
!            !  Line repeated - otherwise could have failure if
!            !  longitudes W are > pi rather than < 0.
             OMEGAE = OMEGAB + DTRAD
             IF (OMEGAE >  TWOPI) OMEGAE = OMEGAE - TWOPI
             OMEGAS = 2. * HLD
! Now that the start-time, end-time and sunset are set in terms of hour
! angle, can set the two hour-angles for (3.2.3).  The simple cases are
! start-to-end-of-timestep, start-to-sunset, sunrise-to-end and sunrise-
! -to-sunset, but two other cases exist and need special treatment.
             IF (OMEGAB <= OMEGAS .OR. OMEGAB <  OMEGAE) THEN
                OMEGA1 = OMEGAB - HLD
              ELSE
                OMEGA1 = - HLD
             ENDIF
             IF (OMEGAE <= OMEGAS) THEN
                OMEGA2 = OMEGAE - HLD
              ELSE
                OMEGA2 = OMEGAS - HLD
             ENDIF
             IF (OMEGAE >  OMEGAB.AND.OMEGAB >  OMEGAS) OMEGA2=OMEGA1
!  Put in an arbitrary marker for the case when the sun does not rise
!  during the timestep (though it is up elsewhere at this latitude).
!  (Cannot set COSZ & LIT within the ELSE ( COSHLD < 1 ) block
!  because 3.2.3 is done outside this block.)
          ENDIF           ! This finishes the ELSE (perpetual day) block
          DIFSIN = SIN(OMEGA2) - SIN(OMEGA1)             ! Begin (3.2.3)
          DIFTIM = OMEGA2 - OMEGA1
          MEAN_OMEGA(J) = (OMEGA1 + OMEGA2)/2.
! Next, deal with the case where the sun sets and then rises again
! within the timestep.  There the integration has actually been done
! backwards over the night, and the resulting negative DIFSIN and DIFTIM
! must be combined with positive values representing the whole of the
! timestep to get the right answer, which combines contributions from
! the two separate daylit periods.  A simple analytic expression for the
! total sun throughout the day is used.  (This could of course be used
! alone at points where the sun rises and then sets within the timestep)
          IF (DIFTIM <  0.) THEN
            DIFSIN = DIFSIN + 2. * SQRT(1.-COSHLD**2)
            DIFTIM = DIFTIM + 2. * HLD
            MEAN_OMEGA(J) = MEAN_OMEGA(J) + Pi
          ENDIF
          IF (MEAN_OMEGA(J) >  Pi) MEAN_OMEGA(J)=MEAN_OMEGA(J)-TWOPI
          IF (DIFTIM == 0.) THEN
! Pick up the arbitrary marker for night points at a partly-lit latitude
             COSZ(J) = 0.
             LIT(J) = 0.
           ELSE
             COSZ(J) = DIFSIN*COSCOS/DIFTIM + SINSIN     ! (3.2.3)
             LIT(J) = DIFTIM / DTRAD
          ENDIF
       ENDIF            ! This finishes the ELSE (perpetual night) block
  100 CONTINUE

      IF (ALLOCATED(sol_bearing)) THEN

        LAT = ASIN(SINLAT)
        DEC = ASIN(SINDEC)

        X = ATAN2( -COS(DEC)*SIN(MEAN_OMEGA),                           &
     &        COS(LAT)*SINDEC - SINLAT*COS(DEC)*COS(MEAN_OMEGA) )

        sol_bearing=RESHAPE( MODULO(X , TWOPI) ,                        &
     &              SHAPE(sol_bearing))

      END IF

      RETURN
      END SUBROUTINE SOLANG
#endif
