#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************
!
!+ Subroutine to set dimensions for the radiation code.
!
! Purpose:
!   To set the dimensions of arrays for use in the radiation code
!   depending on the options chosen.
!
! Method:
!   This routine covers a range of disparate requirements, working by
!   IF-tests.
!
! Current owner of code: James Manners
!
! Description of code:
!   Fortran 90
!
!- ---------------------------------------------------------------------
SUBROUTINE r2_set_rad_dim &
!
(row_length, rows, sin_lat, lon, l_solar, sindec, h_angle, &
 l_geostationary, sat_hgt, sat_lat, sat_lon, &
 i_cloud, i_angular_integration, i_sph_mode, &
 nd_field, n_points, model_levels, cloud_levels, &
 nd_cloud_type, nd_cloud_component, &
 l_extra_top, id_ct, n_rad_layers, nd_column, &
 nd_band, n_band, map_channel, n_channel, nd_channel, &
 nd_viewing_level, n_viewing_level, viewing_level, &
 nd_direction, n_view_direction, view_direction, &
 nd_brdf_basis_fnc, nd_brdf_trunc, &
 nd_profile, nd_flux_profile, nd_radiance_profile, &
 nd_field_flux_diag, nd_field_rad_diag, &
 l_MOSES_ii, l_ctile, nd_tile, nd_point_tile &
)
!
!
!
  IMPLICIT NONE
!
!
! Header files.
#include "cloud_scheme_pcf3z.h"
#include "angular_integration_pcf3z.h"
#include "sph_mode_pcf3z.h"
#include "c_a.h"
!
!
! Dummy arguments
!
  INTEGER, Intent(IN) :: row_length
!   Length of EW rows
  INTEGER, Intent(IN) :: rows
!   Number of NS rows
!
  REAL, Dimension(row_length, rows) :: sin_lat
!   Sines of latitudes of grid-points
  REAL, Dimension(row_length, rows) :: lon
!   Longitudes of grid-points
  LOGICAL, Intent(IN) :: l_solar
!   Flag for solar region of the spectrum
  REAL, Intent(IN) :: sindec
!   Sine of solar declination
  REAL, Intent(IN) :: h_angle(row_length,rows)
!   Hour Angle

!
  INTEGER, Intent(IN) :: i_cloud
!   Option for cloud geomtery
  INTEGER, Intent(IN) :: i_angular_integration
!   Type of angular integration
  INTEGER, Intent(IN) :: i_sph_mode
!   Mode of operation of the spherical harmonic code
  LOGICAL, Intent(IN) :: l_extra_top
!   Flag to include an extra layer in radiation to model the part
!   of the atmosphere above the top of the model
!
  LOGICAL, Intent(IN) :: l_MOSES_ii
!   Flag to use the MOSES II surface scheme
  LOGICAL, Intent(IN) :: l_ctile
!   Flag for coastal tiling
!
  INTEGER, Intent(IN) :: nd_field
!   Horizontal size of input fields
  INTEGER, Intent(IN) :: nd_band
!   Size allocated for spectral bands
  INTEGER, Intent(IN) :: n_band
!   Number of spectral bands used
  INTEGER, Intent(IN) :: n_points
!   Number of grid points to operate on
  INTEGER, Intent(IN) :: model_levels
!   Number of theta levels in the model
  INTEGER, Intent(IN) :: cloud_levels
!   Number of potentially cloudy layers in the model
!
  LOGICAL, Intent(IN) :: l_geostationary
!   Flag for a geostationary orbit
  REAL, Intent(IN) :: sat_hgt
!   Height of the satellite's orbit
  REAL, Intent(IN) :: sat_lon
!   Longitude of the (geostationary) satellite
  REAL, Intent(IN) :: sat_lat
!   Laitude of the (geostationary) satellite
!
  INTEGER, Intent(OUT) :: n_channel
!   Number of channels in output
  INTEGER, Intent(OUT), Dimension(nd_band) :: map_channel
!   Mapping of spectral bands to output channels
!
  INTEGER, Intent(OUT) :: n_rad_layers
!   Number of layers used in radiation
  INTEGER, Intent(OUT) :: nd_cloud_type
!   Size to be allocated for types of cloud
  INTEGER, Intent(OUT) :: nd_cloud_component
!   Size to be allocated for condensed components
  INTEGER, Intent(OUT) :: nd_column
!   Size to be allocated for subcolumns in a grid-box
  INTEGER, Intent(OUT) :: id_ct
!   Topmost layer allowed to contain cloud
  INTEGER, Intent(OUT) :: nd_channel
!   Size to be allocated for satellite channels available from
!   one call to the radiation scheme
  INTEGER, Intent(OUT) :: nd_viewing_level
!   Size to be allocated for viewing levels, where fluxes or radiances
!   will be calculated
  INTEGER, Intent(OUT) :: nd_direction
!   Size to be allocated for viewing directions
  INTEGER, Intent(OUT) :: nd_profile
!   Size to be allocated for points where fields are required
!   for both flux and radiance calculations
  INTEGER, Intent(OUT) :: nd_radiance_profile
!   Size to be allocated for points, where radiances
!   will be calculated
  INTEGER, Intent(OUT) :: nd_flux_profile
!   Size to be allocated for points, where fluxes
!   will be calculated
  INTEGER, Intent(OUT) :: nd_field_flux_diag
!   Size to be allocated for flux diagnostics
  INTEGER, Intent(OUT) :: nd_field_rad_diag
!   Size to be allocated for radiance diagnostics
!
  INTEGER, Intent(OUT) :: nd_brdf_basis_fnc
!   Size to be allocated for basis functions for BRDFs
  INTEGER, Intent(OUT) :: nd_brdf_trunc
!   Size to be allocated for order of truncation fo be applied to
!   BRDFs
!
  INTEGER, Intent(OUT) :: nd_tile
!   Size to be allocated for number of surface tiles within the
!   radiation scheme
  INTEGER, Intent(OUT) :: nd_point_tile
!   Size to be allocated for number of points where surface tiling will
!   be applied
!
  INTEGER, Intent(OUT) :: n_viewing_level
!   Actual number of viewing levels
  REAL, Dimension(model_levels+1)  :: viewing_level
!   Viewing levels (this hard-wired dimension is a temporary measure)
!
  INTEGER, Intent(OUT) :: n_view_direction
!   Actual number of viewing directions
  REAL, Dimension(nd_field, 1, 2)  :: view_direction
!   Viewing directions (this hard-wired dimension is a temporary measure)
!   The second dimension allows just one viewing direction
!
!
!
! Local variables
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: l
!   Loop variable
!
  REAL, Dimension(3) :: r_sat
!   Cartesian coordinates of the satellite
  REAL, Dimension(3) :: r_obs
!   Cartesian coordinates of observed point
  REAL, Dimension(3) :: r_diff
!   Cartesian coordinates of separation of the points
  REAL, Dimension(3) :: r_sun
!    Cartesian coordinates of the sun
  REAL, Dimension(3) :: r_hor
!    Projection of r_diff on local horizontal plane
  REAL, Dimension(3) :: r_hor_s
!    Projection of r_sun on local horizontal plane
  REAL :: mag_rdiff
!   Magnitude of the separation
  REAL :: mag_robs
!   Magnitude of r_obs
  REAL :: mag_rhor
!    Magnitude or r_hor
  REAL :: mag_rhor_s
!    Magnitude of r_hor_s
  REAL :: solar_zen
!    Solar zenith angle
  REAL :: pi_val
!
!
!
!
!
! Set the actual size of arrays in the radiation code:
! for some architectures (such as that of Cray vector
! machines) on odd size is preferred to avoid memory
! bank conflicts.
  nd_profile=2*(n_points/2)+1
!
! Set the number of layers seen in the radiation code.
! This may optionally be 1 greater than the number used
! in the rest of the model to avoid spurious effects
! resulting from the upper boundary (principally in
! stratospheric configurations).
  if (l_extra_top) then
    n_rad_layers=model_levels+1
  else
    n_rad_layers=model_levels
  endif
!
! Standard sizes for cloud arrays.
  nd_cloud_type      = 4
  nd_cloud_component = 4
  SELECT CASE(i_cloud)
    CASE(IP_cloud_column_max)
      nd_column = 3 * cloud_levels + 2
    CASE DEFAULT
      nd_column = 1
  END SELECT
!
!
  IF (   (i_angular_integration == IP_two_stream) .OR. &
       ( (i_angular_integration == IP_spherical_harmonic) .AND. &
         (i_sph_mode == IP_sph_mode_flux) ) ) THEN
!
    SELECT CASE(i_angular_integration)
      CASE(IP_two_stream)
        nd_viewing_level     = 1
        n_viewing_level      = 0
        nd_radiance_profile  = 1
      CASE(IP_spherical_harmonic)
        nd_viewing_level     = model_levels+1
        n_viewing_level      = model_levels+1
        nd_radiance_profile  = nd_profile
    END SELECT
    nd_direction       = 1
    nd_brdf_basis_fnc  = 2
    nd_brdf_trunc      = 1
    nd_field_flux_diag = nd_field
    nd_flux_profile    = nd_profile
    nd_field_rad_diag  = 1
!
    DO i=1, model_levels + 1
      viewing_level(i) = REAL(i-1)
    ENDDO
!
!   Dummy initailization.
    n_view_direction = 1
!
!   Set the mapping of bands in the spectral file onto channels for
!   the diagnostics. When calculating fluxes all bands are normally
!   mapped into a single channel.
    nd_channel            = 1
    n_channel             = 1
    map_channel(1:n_band) = 1
!
  ELSE
!
!   Set up space for only one viewing direction initially.
!
    nd_direction     = 1
    n_view_direction = 1
!
!   Set the mapping of bands in the spectral file onto channels for
!   the diagnostics. We assume for now that each band in the spectral
!   file corresponds to a particular channel.
    nd_channel = n_band
    n_channel  = n_band
    DO i=1, n_band
      map_channel(i)=i
    ENDDO
!
    SELECT CASE(l_geostationary)
!
      CASE(.TRUE.)
!
! Define the Cartesian position of the satellite.

        r_sat(1) = ( Earth_radius + sat_hgt ) * &
                     COS(sat_lat) * COS(sat_lon)
        r_sat(2) = ( Earth_radius + sat_hgt ) * &
                     COS(sat_lat) * SIN(sat_lon)
        r_sat(3) = ( Earth_radius + sat_hgt ) * &
                     SIN(sat_lat)

!
        DO j=1, rows
          DO i=1, row_length
!
            l = i+(j-1)*row_length

!
! Calculate zenith angle
!

! Define point in satelite footprint the satelite is lookin at:

            r_obs(1) = Earth_radius * &
                         COS(ASIN(sin_lat(i, j))) * COS(lon(i, j))
            r_obs(2) = Earth_radius * &
                         COS(ASIN(sin_lat(i, j))) * SIN(lon(i, j))
            r_obs(3) = Earth_radius * &
                         SIN(ASIN(sin_lat(i, j)))
!
            r_diff(1) = r_sat(1) - r_obs(1)
            r_diff(2) = r_sat(2) - r_obs(2)
            r_diff(3) = r_sat(3) - r_obs(3)
!
            mag_rdiff = SQRT( r_diff(1) * r_diff(1) + &
                              r_diff(2) * r_diff(2) + &
                              r_diff(3) * r_diff(3) )
!
            mag_robs  = SQRT( r_obs(1) * r_obs(1) + &
                              r_obs(2) * r_obs(2) + &
                              r_obs(3) * r_obs(3) )

!
            view_direction(l, 1, 1) = ( r_obs(1) * r_diff(1) + &
                                        r_obs(2) * r_diff(2) + &
                                        r_obs(3) * r_diff(3) ) / &
                                      ( Earth_Radius * mag_rdiff )

!            pi_val=4.0*DATAN(1.0)
!            view_direction(l, 1, 1)=pi_val/3.40
!
! Calculate Azimuthal angle
!
            IF (l_solar) THEN

! Define the normalised vecto pointing to the sun

              r_sun(1)=COS(sindec)*COS(h_angle(i,j))
              r_sun(2)=COS(sindec)*SIN(h_angle(i,j))
              r_sun(3)=SIN(sindec)

! Projection of r_diff on local horizontal plane

              r_hor(1)=r_diff(1)-view_direction(l,1,1)*r_diff(1)
              r_hor(2)=r_diff(2)-view_direction(l,1,1)*r_diff(2)
              r_hor(3)=r_diff(3)-view_direction(l,1,1)*r_diff(3)

!
! Projection of unity vector r_sun on horizontal plane
!

              solar_zen=(r_obs(1)*r_sun(1)+r_obs(2)*r_sun(2)  &
                       +r_obs(3)*r_sun(3))/mag_robs

              r_hor_s(1)=r_sun(1)-solar_zen*r_sun(1)
              r_hor_s(2)=r_sun(2)-solar_zen*r_sun(2)
              r_hor_s(3)=r_sun(3)-solar_zen*r_sun(3)

!
! Calculate the azimuth angle and store in view_direction(l,1,2)
!

              mag_rhor  = SQRT( r_hor(1) * r_hor(1) + &
                                r_hor(2) * r_hor(2) + &
                                r_hor(3) * r_hor(3) )

              mag_rhor_s  = SQRT( r_hor_s(1) * r_hor_s(1) + &
                                  r_hor_s(2) * r_hor_s(2) + &
                                  r_hor_s(3) * r_hor_s(3) )

              view_direction(l,1,2)=(r_hor(1)*r_hor_s(1) &
                                    +r_hor(2)*r_hor_s(2) &
                                    +r_hor(3)*r_hor_s(3)) &
                                    /mag_rhor/mag_rhor_s

!              view_direction(l, 1, 2) = 0.0

!
            ELSE
!
!             For infra-red simulations we set the azimuth to 0.0.
              view_direction(l, 1, 2) = 0.0
!
            ENDIF
!
          ENDDO
        ENDDO
!
      CASE(.FALSE.)
!
! DEPENDS ON: ereport
        CALL Ereport("rddim", 21, &
          "Only geostationary satellites are available so far")
!
    END SELECT
!
!   There should be only one viewing level: that at the top of the
!   atmosphere.
    nd_viewing_level  = 1
    n_viewing_level  = 1
    viewing_level(1) = 0.0
!
    nd_field_flux_diag  = 1
    nd_flux_profile     = 1
    nd_radiance_profile = nd_profile
    nd_brdf_basis_fnc   = 2
    nd_brdf_trunc       = 0
    nd_field_flux_diag  = 1
    nd_field_rad_diag   = nd_field
!
  ENDIF
!
! Remaining common sizes.
  nd_point_tile = n_points
  IF (l_MOSES_ii .AND. l_ctile) THEN
    nd_tile     = 3
  ELSE
    nd_tile     = 2
  ENDIF
  id_ct         = n_rad_layers+1-cloud_levels
!
!
!
END SUBROUTINE r2_set_rad_dim
#endif
