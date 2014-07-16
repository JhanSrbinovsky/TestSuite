#if defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Longwave interface to the radiance code.
!
! Purpose:
!   This routine prepares the call to the ES radiance
!   scheme in the longwave.
!
! Method:
!   Principally, this routine transfers arrays into the correct formats.
!
! Current owner of code: J.-C. Thelen
!
! History:
!  Ver   Date      Comment
!  6.2   25/01/06  Code included into UM-build 6.2 (Jean-Claude Thelen)
!  6.2   25/01/06  Pass l_inhom_cloud and inhom_cloud to radiance_calc.
!                                                   (James Manners)
!
! Description of code:
!   FORTRAN 90
!
!- ---------------------------------------------------------------------
!+ Subroutine to set surface fields.
!
! Purpose:
!   The albedos and emissivity of the surface are set.
!
! Method:
!   Straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             08-03-06                Original code
!                                               (J.-C. Thelen)
!
! Description of code:
!   FORTRAN 90
!
!- ---------------------------------------------------------------------
      subroutine r2_set_surface_field_lw(ierr                           &
     &  , n_points, list, n_band, ls_brdf_trunc                         &
     &  , l_ctile, flandg                                               &
     &  , n_brdf_basis_fnc, f_brdf, rho_alb                             &
     &  , land, ice_fraction, tstar_sea, tstar_solid                    &
     &  , l_rad_tile, n_point_tile, n_tile, list_tile, rho_alb_tile     &
     &  , frac_tile, t_tile                                             &
     &  , list_tile_outer, index_tile                                   &
     &  , nd_field, nd_profile, nd_band, nd_brdf_basis_fnc              &
     &  , nd_brdf_trunc, nd_point_tile, nd_tile                         &
     &  )
!
!
!
      use tileid3z
!
      implicit none
!
!
!     Comdecks included
#include "c_0_dg_c.h"
#include "def_std_io_icf3z.h"
#include "surface_spec_pcf3z.h"
#include "error_pcf3z.h"
!
!     Dummy variables:
!
!     Dimensions of arrays:
      integer                                                           &
                !, intent(in)
     &    nd_field                                                      &
!           Allocated size of fields of data
     &  , nd_profile                                                    &
!           Maximum number of atmospheric profiles
     &  , nd_band                                                       &
!           Maximum number of spectral bands
     &  , nd_brdf_basis_fnc                                             &
!           Maximum number of BRDF basis functions
     &  , nd_brdf_trunc                                                 &
!           Maximum order of truncation of BRDF
     &  , nd_point_tile                                                 &
!           Size allocated for points where the surface is tiled
     &  , nd_tile
!           Size allocated for surface tiles
!
      integer                                                           &
                !, intent(out)
     &    ierr
!           Error flag
!
!     Actual sizes used:
      integer                                                           &
                !, intent(in)
     &    n_points                                                      &
!           Number of atmospheric points
     &  , n_band
!           Number of spectral bands
!
!     Points to be treated:
      integer                                                           &
                !, intent(in)
     &    list(nd_field)
!           List of sunlit points
!
      logical                                                           &
                !, intent(in)
     &    land(nd_field)
!           Land flag
      LOGICAL, Intent(IN) :: l_ctile
!           Flag for coastal tiling
      real                                                              &
                !, intent(in)
     &    flandg(nd_field)
!           land fraction in grid box
      real                                                              &
                !, intent(in)
     &    ice_fraction(nd_field)                                        &
!           Ice fractions in oceanic grid-boxes
     &  , tstar_sea(nd_field)                                           &
!           Temperatures of the sea surface
     &  , tstar_solid(nd_field)
!           Temperatures of the solid surface
!
!     Physical properties of surfaces:
      integer                                                           &
                !, intent(in)
     &    ls_brdf_trunc
!           Order of truncation applied to BRDFs
!
!     Surface properties set.
      integer                                                           &
                !, intent(out)
     &    n_brdf_basis_fnc
!           Number of basis functions for BRDFs
      real                                                              &
                !, intent(out)
     &    f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                  &
     &      , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)                     &
!           Basis functions for the surface
     &  , rho_alb(nd_profile, nd_brdf_basis_fnc, nd_band)
!           Weights of the surface BRDF basis functions
!
!     Arrays related to tiling of the surface
      logical                                                           &
                !, intent(in)
     &    l_rad_tile
!           Local to allow tiling options
      integer                                                           &
                !, intent(out)
     &    n_point_tile                                                  &
!           Number of points to tile
     &  , n_tile                                                        &
!           Number of tiles used
     &  , list_tile(nd_point_tile)                                      &
!           List of points with surface tiling
     &  , list_tile_outer(nd_point_tile)                                &
!           List of points with surface tiling in the full list
     &  , index_tile(npd_tile_type)
!           The indexing number of tiles of the given type
      real                                                              &
                !, intent(out)
     &    rho_alb_tile(nd_point_tile, nd_brdf_basis_fnc                 &
     &      , nd_tile, nd_band)                                         &
!           Weights for the basis functions of the BRDFs
!           at the tiled points
     &  , t_tile(nd_point_tile, nd_tile)                                &
!           Local surface temperatures on individual tiles
     &  , frac_tile(nd_point_tile, nd_tile)
!           Fractions of each tiled grid-point occupied by tiles
!           of the appropriate type
!
!
!     Local variables.
!
      integer                                                           &
     &    i                                                             &
!           Loop variable
     &  , l                                                             &
!           Loop variable
     &  , ll
!           Loop variable
!
!
!
!     Zero the surface albedos.
      do i=1, n_band
        rho_alb(1:n_points, ip_surf_alb_diff, i) = 0.0
        rho_alb(1:n_points, ip_surf_alb_dir, i) = 0.0
      enddo
!
!     Set the surface basis functions for a Lambertian surface.
      n_brdf_basis_fnc=1
!     By defining F_{1,0,0,0} to be 4, RHO_ALB beomes equal to the
!     diffuse albedo.
      f_brdf(1, 0, 0, 0)=4.0
      if (ls_brdf_trunc /= 0) then
        write(iu_err, '(/a)')                                           &
     &    'Error: The order of surface truncation is too high.'
        ierr=i_err_fatal
        return
      endif
!
!
      if (l_rad_tile) then
!
!       Set up the surface tiling variables. There are multiple levels
!       of indexing. Over all points in the domain only those in the
!       array list require radiative calculations and of these points
!       only those in the array list_file require tiling, this array
!       being indexed over points where radiative calculations are to
!       be done. list_tile_outer is indexed over the tiled points and
!       gives the index in the whole domain.
!
        if (l_ctile) then
!
!         With coastal tiling we can have land, open sea or sea ice
!         in the grid-box. Coastal tiling is not possible without
!         MOSES-II.
!
          n_tile=3
          index_tile(ip_ocean_tile)=1
          index_tile(ip_seaice_tile)=2
          index_tile(ip_land_tile)=3
          n_point_tile=0
          do ll=1, n_points
            l=list(ll)
            if ( (flandg(l) < 1.0) .AND.                                &
     &           ( (flandg(l) > 0.0) .OR.                               &
     &             ( (ice_fraction(l) > 0.0) .AND.                      &
     &               (ice_fraction(l) < 1.0) ) ) ) then
              n_point_tile=n_point_tile+1
              list_tile(n_point_tile)=ll
              list_tile_outer(n_point_tile)=l
!             Assign tiled fractions consistent with the indices above.
              frac_tile(n_point_tile, 1)                                &
     &          =(1.0-flandg(l))*(1.0-ice_fraction(l))
              frac_tile(n_point_tile, 2)                                &
     &          =(1.0-flandg(l))*ice_fraction(l)
              frac_tile(n_point_tile, 3)                                &
     &          =flandg(l)
            endif
          enddo
!
        else
!
!         Without coastal tiling we have only open sea or sea ice
!         forming the coastal tiling. The setting of the land albedos
!         depends on whether MOSES-II is on.
!
          n_tile=2
          index_tile(ip_ocean_tile)=1
          index_tile(ip_seaice_tile)=2
          n_point_tile=0
          do ll=1, n_points
            l=list(ll)
            if ( (.not.land(l)).and.                                    &
     &           (ice_fraction(l) >  0.0).and.                          &
     &           (ice_fraction(l) <  1.0) ) then
              n_point_tile=n_point_tile+1
              list_tile(n_point_tile)=ll
              list_tile_outer(n_point_tile)=l
              frac_tile(n_point_tile, 1)                                &
     &          =(1.0-ice_fraction(l))
              frac_tile(n_point_tile, 2)                                &
     &          =ice_fraction(l)
            endif
          enddo
!
        endif
!
!       Note that if non-zero surface emissivities are used the
!       following block will need to be revised along the lines of
!       the SW.
!
        rho_alb_tile(1:n_point_tile, ip_surf_alb_dir                    &
     &    , 1:n_tile, 1:n_band)  = 0.0
        rho_alb_tile(1:n_point_tile, ip_surf_alb_diff                   &
     &    , 1:n_tile, 1:n_band) = 0.0
!
!       Tiled temperatures are required to handle emission from the
!       surface. Ensure that the indexing of the second subscript is
!       consistent with the assignment of index_tile above.
        t_tile(1:n_point_tile, 1)                                       &
     &    = tstar_sea(list_tile_outer(1:n_point_tile))
        t_tile(1:n_point_tile, 2)                                       &
     &    = tstar_solid(list_tile_outer(1:n_point_tile))
        if (l_ctile) then
          t_tile(1:n_point_tile, 3)                                     &
     &      = tstar_solid(list_tile_outer(1:n_point_tile))
        endif
!
      endif
!
!
!
      return
      END SUBROUTINE r2_set_surface_field_lw
!+ Subroutine to set algorithmic options.
!
! Purpose:
!   Algorithmic options and array sizes to be set interactively
!   are determined.
!
! Method:
!   Straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             08-03-06                Original code
!                                               (J.-C. Thelen)
!
! Description of code:
!   FORTRAN 90
!
!- ---------------------------------------------------------------------
#endif
