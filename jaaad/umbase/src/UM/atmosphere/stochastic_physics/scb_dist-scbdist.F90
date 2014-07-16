#if defined(A35_1A)
      REAL FUNCTION SCB_DIST(LAT1,LAT2,LONG1,LONG2)
!-----------------------------------------------------------------------
! THE DISTANCE BETWEEN ANY TWO POINTS ON A SPHERE IS CALCULATED FROM
!   COS(D) = SIN(A1)SIN(A2) + COS(A1)COS(A2)COS(B1-B2)
! WHERE A IS THE LATITUDE, B IS THE LONGITUDE AND D IS THE ANGLE BETWEEN
! THE TWO POINTS.
! THIS WILL BE USED TO WORK OUT THE DISTANCE BETWEEN TWO POINTS.
!-----------------------------------------------------------------------
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  6.2  28/02/05   Function to calculate great circle distances
!                  needed for STPH_SCV. Original code by M.Gray
!                  Adapted and included by A.Arribas
!
#include "c_a.h"
!
      REAL, INTENT(InOut) :: LAT1,LAT2,LONG1,LONG2
      REAL :: COSD  ! Local Variable
                    ! It stores the distance between two points
!
      IF(LAT1 == LAT2 .AND. LONG1 == LONG2)THEN
        SCB_DIST=0.0
      ELSE
        COSD = (SIN(LAT1)*SIN(LAT2)) +                                  &
     &         (COS(LAT1)*COS(LAT2)*COS(LONG1-LONG2))
        IF( COSD >  1.0 )COSD=1.0
        IF( COSD <  -1.0)COSD=-1.0
        SCB_DIST=ABS( ACOS(COSD)*Earth_Radius )
      ENDIF
!
      END FUNCTION SCB_DIST
#endif
