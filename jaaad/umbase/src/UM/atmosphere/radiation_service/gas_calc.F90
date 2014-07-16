#if defined(A70_1A) || defined(A70_1B) \
 || defined(A70_1C) || defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL Subroutine GAS_CALC ----------------------------------------------
!LL
!LL Purpose :
!LL   Calculates the trace gas mixing ratio (or weighting factor for
!LL aerosol forcing fields.  Rates of increase (yearly compound factors)
!LL can be supplied, or spot values (which will be linearly
!LL interpolated) or a mixture of these.  It is designed so it can be
!LL called each time step, but when rates of increase are being used,
!LL values are in fact only updated at New Year.
!LL The rules are:
!LL   If rates exist (i.e. are positive) for the first & current years
!LL then all concentrations are ignored, except for the initial value.
!LL   If there is a positive rate for the current year but not for the
!LL start, the current rate & most recent concentration are used.
!LL   If rates do not exist for the current year then the concentration
!LL is calculated by linear interpolation between the concentrations at
!LL the given years.
!LL   The mixing ratios calculated after the last given year use the
!LL rate for the final given year.
!LL   The 360-day year is assumed.
!LL CARE should be taken if solitary rates are specified, as this can
!LL result in discontinuities in the concentration time profile at
!LL the next given year without a corresponding given rate.
!LL
!LL Authors : Andrew Brady, Tim Johns, William Ingram
!LL
!LL Version for : Cray YMP
!LL
!LL  Model
!LL version  Date
!LL   4.2  19/11/96       Author: William Ingram, reviewer Cath Senior.
!LL   4.4  22/9/97         Correct code for the case when one changes
!LL    from linear interpolation to a rate & then changes the rate.  WJI
! 6.2  21/02/06     Introduced #if defined statements for
!                   version control of radiation code.
!                                  (J.-C.Thelen)
!LL
!LL
!LLEND -----------------------------------------------------------------
!*L Arguments

      SUBROUTINE GAS_CALC(GAS_NOW                                       &
     &                   ,GAS_INDEX_MAX                                 &
     &                   ,GAS_YEAR                                      &
     &                   ,GAS_CONC                                      &
     &                   ,GAS_RATE                                      &
     &                   ,MAX_SCENARIO_PTS                              &
     &                   ,ICODE)

      IMPLICIT NONE

! Maximum size of common block arrays
#include "cmaxsize.h"
! Submodel parameters for array sizes
#include "csubmodl.h"

      REAL          GAS_NOW      !OUT Gas concentration at time step
      INTEGER       GAS_INDEX_MAX                                       &
                                 !IN
     &            , MAX_SCENARIO_PTS ! IN
      INTEGER       GAS_YEAR(MAX_SCENARIO_PTS)   !IN
      REAL          GAS_CONC(MAX_SCENARIO_PTS)   !IN
      REAL          GAS_RATE(MAX_SCENARIO_PTS)   !IN
      INTEGER       ICODE        !OUT Return code: successful=0
!*

! Common blocks

! Model time
#include "ctime.h"

!     Local variables

      INTEGER       INDEX       ! to subscript gas concs for NOW_TIME
      INTEGER       I           ! Loop over indices
      INTEGER       YEAR_IN_SECS! Year length in seconds
      INTEGER       NOW_TIME_DAY, NOW_TIME_SEC
!                               ! Time now in days/secs from time zero
      INTEGER       GAS_YR_DAY1,  GAS_YR_SEC1
!                               ! Time in days/secs of current GAS_YEAR
      INTEGER       TIME1
!                               ! The same converted to seconds
      INTEGER       GAS_YR_DAY2,  GAS_YR_SEC2
!                               ! Time in days/secs of next GAS_YEAR

!     Check that GASCNST namelist is defined for this year
      IF ( I_YEAR  <   GAS_YEAR(1) ) THEN
        ICODE = 8325
        WRITE (6, *) 'GAS_CALC: no gas data for this year' 
        RETURN
      ENDIF

!     Loop over I to find correct index for current NOW_TIME
      INDEX = 0
      DO I=1, GAS_INDEX_MAX
        IF ( I_YEAR  >=  GAS_YEAR(I) ) INDEX = INDEX+1
      ENDDO

!     Calculate time now in seconds
! DEPENDS ON: time2sec
      CALL TIME2SEC (I_YEAR, I_MONTH, I_DAY, I_HOUR, I_MINUTE, I_SECOND,&
     &              0, 0, NOW_TIME_DAY, NOW_TIME_SEC, .TRUE.)

!     If gas rate at current year is non zero calculate new GAS_NOW
!     by considering compound increases of GAS_RATE(1:INDEX)
      IF ( GAS_RATE(INDEX)  >   0. ) THEN
        YEAR_IN_SECS = 360 * 86400
! DEPENDS ON: time2sec
        CALL TIME2SEC (GAS_YEAR(INDEX), 1, 1, 0, 0, 0,                  &
     &                0, 0, GAS_YR_DAY1, GAS_YR_SEC1, .TRUE.)
        GAS_NOW = GAS_CONC(1)
        DO I=1, INDEX-1
          IF ( GAS_RATE(I)  <   0. ) THEN
             GAS_NOW = GAS_CONC(I+1)
           ELSE
             GAS_NOW = GAS_NOW *                                        &
     &            ( GAS_RATE(I) ** REAL(GAS_YEAR(I+1)-GAS_YEAR(I)) )
          ENDIF
        ENDDO
!       GAS_NOW now holds the concentration in year INDEX - need only
!       update it to the current year.
        GAS_NOW=GAS_NOW*(GAS_RATE(INDEX)**                              &
     &    REAL(((NOW_TIME_DAY-GAS_YR_DAY1)*86400+                       &
     &          NOW_TIME_SEC-GAS_YR_SEC1)/YEAR_IN_SECS))

!     Otherwise calculate by linear interpolation between respective
!     GAS concentrations of given years.
      ELSE
! DEPENDS ON: time2sec
        CALL TIME2SEC (GAS_YEAR(INDEX), 1, 1, 0, 0, 0,                  &
     &                0, 0, GAS_YR_DAY1, GAS_YR_SEC1, .TRUE.)
! DEPENDS ON: time2sec
        CALL TIME2SEC (GAS_YEAR(INDEX+1), 1, 1, 0, 0, 0,                &
     &                0, 0, GAS_YR_DAY2, GAS_YR_SEC2, .TRUE.)
        TIME1   = GAS_YR_DAY1*86400 + GAS_YR_SEC1
        GAS_NOW = GAS_CONC(INDEX) +                                     &
     &          ( GAS_CONC(INDEX+1) - GAS_CONC(INDEX) )                 &
     & * REAL ( NOW_TIME_DAY*86400 + NOW_TIME_SEC - TIME1 )             &
     &      / REAL ( GAS_YR_DAY2*86400 + GAS_YR_SEC2 - TIME1 )
      ENDIF

      RETURN
      END SUBROUTINE GAS_CALC
#endif
