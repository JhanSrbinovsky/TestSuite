#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
FUNCTION in_footprint(i_call, l_sw, &
  l_geostationary, min_view_lon, max_view_lon, &
  min_view_lat, max_view_lat, &
  Sin_Pt_lat, Pt_lon, time_0) &
!
RESULT (l_in)
!
! Purpose: Routine to determine whether a grid-point is within the
!          footprint seen by a satellite.
!
! Method:  This routine is subject to considerable alteration and
!          will need rewriting for each application. A satellite is
!          identified only by the number of the call.
!
! History:
! Version   Date     Comment
! ---      -----     -------
! 6.2      13/02/06  Original code included into UM Build
!                              (J.-C. Thelen)
!
! Code description:
!   Language: Fortran 90
!
!
! Declarations
!
  IMPLICIT NONE
!
!
! Header files.
#include "c_pi.h"
!
!
  INTEGER, Intent(IN) :: i_call
!   Number of call to radiation code
  LOGICAL, Intent(IN) :: l_sw
!   Flag for Shortwave region
  REAL, Intent(IN) :: Pt_lon
!   Longitude of grid-point
  REAL, Intent(IN) :: Sin_Pt_lat
!   Sine of latitude of grid-point
!   (available directly from calling routine)
  REAL, Intent(IN) :: time_0
!   Time at beginning of current timestep
!
  LOGICAL, Intent(IN) :: l_geostationary
!   Flag for geostationary satellite
  REAL, Intent(IN) :: min_view_lon
!   Minimum viewing longitude
  REAL, Intent(IN) :: max_view_lon
!   Maximum viewing longitude
  REAL, Intent(IN) :: min_view_lat
!   Minimum viewing latitude
  REAL, Intent(IN) :: max_view_lat
!   Maximum viewing latitude
!
  LOGICAL :: l_in
!   Returned value of the function: true if the point is
!   in the footprint
!
!
!
! Local variables
  INTEGER :: Error_code
!   Error code passed to termintaing routine in case of fatal errors
  REAL :: Sat_lon
!   Longitude of viewing satellite
  REAL :: Sat_lat
!   Latitude of viewing satellite
!
  REAL, Parameter :: period =  11000.0
!   Time at beginning of current timestep
!
!
!
! SW calls first:
  SELECT CASE(l_sw)
!
!
    CASE(.TRUE.)
!
      SELECT CASE(i_call)
!
        CASE(1)
!         This should not be called by the main call.
          Error_code = 1
! DEPENDS ON: ereport
          CALL Ereport('in_footprint', Error_code, &
            'Foot-printing in the main call is not permitted.')
        CASE(2)
!
          IF (l_geostationary) THEN
!           Simply test the preset window.
!           Latitudes first:
            l_in = (Sin_pt_lat >= SIN(min_view_lat)) .AND. &
                   (Sin_pt_lat <= SIN(max_view_lat))
!           Longitudes, allowing for different conventions
            l_in = l_in .AND. &
                 ( ( (Pt_lon >= min_view_lon) .AND. &
                     (Pt_lon <= max_view_lon) ) .OR. &
                   ( (Pt_lon - 2.0 * Pi >= min_view_lon) .AND. &
                     (Pt_lon - 2.0 * Pi <= max_view_lon) ) )
          ELSE
!           This code is purely here for illustration.
            Sat_lon = 2.0 * Pi * REAL( time_0 / period - &
                        REAL(INT(time_0 / period) ) )
            Sat_lat = Pi * SIN (Sat_lon)
            IF ( (ABS(Pt_lon - Sat_lon) < 0.2) .AND. &
                 (ABS(Sin_Pt_lat - SIN(Sat_lat)) < 0.2) ) THEN
              l_in = .TRUE.
            ELSE
              l_in = .FALSE.
            ENDIF
          ENDIF
!
        CASE DEFAULT
!         Code is not available for this call.
          Error_code=10
! DEPENDS ON: ereport
          CALL Ereport('in_footprint', Error_code, &
            'No footprinting is available for this call.')
!
      END SELECT
!
!
    CASE(.FALSE.)
!
      SELECT CASE(i_call)
!
        CASE(1)
!         This should not be called by the main call.
          Error_code = 1
! DEPENDS ON: ereport
          CALL Ereport('in_footprint', Error_code, &
            'Foot-printing in the main call is not permitted.')
        CASE(2)
!
          IF (l_geostationary) THEN
!           Simply test the preset window.
!           Latitudes first:
            l_in = (Sin_pt_lat >= SIN(min_view_lat)) .AND. &
                   (Sin_pt_lat <= SIN(max_view_lat))
!           Longitudes, allowing for different conventions
            l_in = l_in .AND. &
                 ( ( (Pt_lon >= min_view_lon) .AND. &
                     (Pt_lon <= max_view_lon) ) .OR. &
                   ( (Pt_lon - 2.0 * Pi >= min_view_lon) .AND. &
                     (Pt_lon - 2.0 * Pi <= max_view_lon) ) )
          ELSE
!           This code is purely here for illustration.
            Sat_lon = 2.0 * Pi * REAL( time_0 / period - &
                        REAL(INT(time_0 / period) ) )
            Sat_lat = Pi * SIN (Sat_lon)
            IF ( (ABS(Pt_lon - Sat_lon) < 0.2) .AND. &
                 (ABS(Sin_Pt_lat - SIN(Sat_lat)) < 0.2) ) THEN
              l_in = .TRUE.
            ELSE
              l_in = .FALSE.
            ENDIF
          ENDIF
!
        CASE DEFAULT
!         Code is not available for this call.
          Error_code=10
! DEPENDS ON: ereport
          CALL Ereport('in_footprint', Error_code, &
            'No footprinting is available for this call.')
!
      END SELECT
!
!
    END SELECT
!
!
!
END FUNCTION in_footprint
#endif
