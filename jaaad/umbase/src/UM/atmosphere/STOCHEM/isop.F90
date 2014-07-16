#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE ISOP(month,isopre)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Calculates multiplier for isoprene emissions
!-                         at particular time.
!-
!-   Inputs  : MONTH,ltdat
!-   Outputs : ISOPRE
!-   Controls:
!-
!
! Current Code Owner: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.5    10/03/95  Created.  W.J. Collins
!  5.5    13/02/04  Vectorised code. K. Ketelsen.
!  6.1    22/10/04  Removed SECS call from inside loop. M.G. Sanderson
!
!-
!VVV  V2.6  ISOP 14/IX/2000  - Instead of ISOPRE
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER,                      INTENT(IN) :: month
      REAL, DIMENSION(nlnpe,nlpe), INTENT(OUT) :: isopre

      INTEGER :: i
      INTEGER :: j
      INTEGER :: hour
      REAL, DIMENSION(0:86400) :: time    !kk size Should be sufficient
      REAL :: xlat
      REAL :: xlong
      REAL :: daylen
      REAL :: thet
      REAL :: ss

!kk   Compute hour outside of loops
! Take time to be 15th of the month
! DEPENDS ON: secs
      ss = SECS(15.0,month,1)
      DO hour = 0, INT(daysec/stochem_advection_step)-1
!       time(hour) = SECS(15.0,month,1) + hour*stochem_advection_step
        time(hour) = ss + hour*stochem_advection_step
      END DO

      DO j=1,nlpe
        DO i=1,nlnpe
          xlat = (lat(j+ltdat-1-1)+lat(j+ltdat-1)) / 2.0 - 90.0
          IF (i+lndat-1 < nlong) THEN
            xlong = (long(i+lndat-1)+long(i+lndat-1+1)) / 2.0
          ELSE
            xlong = (long(i+lndat-1)+long(MOD(i+lndat-1,nlong)+1)+360.0)&
     &        / 2.0
          END IF
          daylen = 0.0
          DO hour = 0, INT(daysec / stochem_advection_step) - 1
! DEPENDS ON: zen
            thet = ZEN(time(hour),xlat,xlong)
            IF (COS(thet) > 0.0) daylen = daylen + COS(thet)
          END DO
          IF(daylen == 0.0) daylen = 0.1
          ISOPRE(I,J) = INT(daysec/stochem_advection_step) / daylen
        END DO
      END DO

      END SUBROUTINE ISOP
#endif
