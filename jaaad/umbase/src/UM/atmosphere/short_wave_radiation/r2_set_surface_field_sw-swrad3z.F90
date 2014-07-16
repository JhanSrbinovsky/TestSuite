#if defined(A01_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set surface fields.
!
! Purpose:
!   The albedos and emissivity of the surface are set.
!
! Method:
!   Straightforward. Though the arrays passed to the code may depend
!   on the spectral band, the input arrays have no spectral dependence.
!   Note that BRDFs are treated as Lambertian here.
!
! current owner of code: James Manners
!
!- ---------------------------------------------------------------------
      subroutine r2_set_surface_field_sw(ierr                           &
     &  , n_band, ls_brdf_trunc                                         &
     &  , nlit, list                                                    &
     &  , l_MOSES_II, l_ctile                                           &
     &  , land, land0p5, open_sea_albedo                                &
     &  , land_alb, sice_alb                                            &
     &  , flandg, ice_fraction                                          &
     &  , land_albedo, weight_690nm                                     &
     &  , land0p5_g, flandg_g                                           &
     &  , n_brdf_basis_fnc, f_brdf, rho_alb                             &
     &  , l_rad_tile, n_point_tile, n_tile, list_tile, rho_alb_tile     &
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
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
#include "surface_spec_pcf3z.h"
!
!     Dummy variables:
!
!     Dimensions of arrays:
      integer                                                           &
                !, intent(in)
     &    nd_field                                                      &
!           Size of input fields
     &  , nd_profile                                                    &
!           Maximum number of atmospheric profiles
     &  , nd_band                                                       &
!           Maximum number of spectral bands
     &  , nd_brdf_basis_fnc                                             &
!           Maximum number of BRDF basis functions
     &  , nd_brdf_trunc                                                 &
!           Maximum order of BRDF terms
     &  , nd_point_tile                                                 &
!           Size allocated for points to be tiled
     &  , nd_tile
!           Size allocated for number of tiles
!
      integer                                                           &
                !, intent(out)
     &    ierr
!           Error flag
!
!     Actual sizes used:
      integer                                                           &
                !, intent(in)
     &    n_band
!           Number of spectral bands
!
!     Lit points:
      integer                                                           &
                !, intent(in)
     &    nlit                                                          &
!           Number of lit points
     &  , list(nd_field)
!           List of sunlit points
!
!     Surface options
      LOGICAL, Intent(IN) :: l_MOSES_II
!       Surface scheme is MOSES-II
      LOGICAL, Intent(IN) :: l_ctile
!       Coastal tiling is used
!
!     Physical properties of surfaces:
      integer                                                           &
                !, intent(in)
     &    ls_brdf_trunc
!           Order of truncation applied to BRDFs
      LOGICAL, Intent(IN) :: land(nd_field)
!           Land mask
      LOGICAL, Intent(IN) :: land0p5(nd_field)
!           Land mask, true where the land fraction exceeds 0.5
      real                                                              &
                !, intent(in)
     &    open_sea_albedo(nd_field, 2)                                  &
!           Diffuse albedo field
     &  , flandg(nd_field)                                              &
!           Fraction of land in a grid-box
     &  , land_alb(nd_field)                                            &
!           Land albedos
     &  , sice_alb(nd_field)                                            &
!           Sea-ice albedos
     &  , land_albedo(nd_field, 4)                                      &
!           Land surface albedo fields for MOSES-II
     &  , weight_690nm(nd_band)                                         &
!           Weights for each band for region below 690 nm
     &  , ice_fraction(nd_field)
!           Fraction of sea ice
!
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
!     Gathered surface fields
      LOGICAL, Intent(OUT), Dimension(nd_profile) :: land0p5_g
!       Gathered land mask: .TRUE. if land fraction > 0.5
      REAL, Intent(OUT), Dimension(nd_profile) :: flandg_g
!       Gathered land fraction
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
     &      , nd_tile, nd_band)
!           Weights for the basis functions of the BRDFs
!           at the tiled points
!
!
!     Local variables.
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
      do l=1, nlit
        land0p5_g(l)=land0p5(list(l))
        flandg_g(l)=flandg(list(l))
      enddo
!
!
!     Define weightings for the basis functions of the surface
!     BRDFs: in effect these are surafce albedos. Each grid-box
!     may contain land, sea or sea-ice, which are treated by tiles
!     within the radiation scheme (which are used more generally
!     then simply when the coastal tiling scheme is on).
!
!     Note: Without coastal tiling the land albedo contains the
!     albedo of the solid surface. MOSES II allows the addition
!     of a spectral dependence, split between the VIS and NIR.
!     If MOMSES-II is selected, coastal tiling may be enabled. In
!     this case, the land alboed refers only to the land surface
!     and there is a separate sea-ice albedo.
!

      do i=1, n_band
        do l=1, nlit
!
!
!         Oceanic surface.
          if (flandg(list(l)) < 1.0) then
            rho_alb(l, ip_surf_alb_diff, i)                             &
     &        =sice_alb(list(l))*ice_fraction(list(l))                  &
     &        +open_sea_albedo(list(l), 2)                              &
     &        *(1.0-ice_fraction(list(l)))
            rho_alb(l, ip_surf_alb_dir, i)                              &
     &        =sice_alb(list(l))*ice_fraction(list(l))                  &
     &        +open_sea_albedo(list(l), 1)                              &
     &        *(1.0-ice_fraction(list(l)))
          else
            rho_alb(l, ip_surf_alb_diff, i) = 0.0
            rho_alb(l, ip_surf_alb_dir, i) = 0.0
          endif
!
!         Add contributions from the land.
          if (flandg(list(l)) > 0.0) then
            if ( l_MOSES_II ) then
!
              if (l_ctile) then
!
                rho_alb(l, ip_surf_alb_diff, i)                         &
     &            =(1.0-flandg_g(l))*rho_alb(l, ip_surf_alb_diff, i)    &
     &            +flandg_g(l) * (weight_690nm(i)                       &
     &            *land_albedo(list(l), 2)                              &
     &            +(1.0-weight_690nm(i))*land_albedo(list(l), 4))
                rho_alb(l, ip_surf_alb_dir, i)                          &
     &            =(1.0-flandg_g(l))*rho_alb(l, ip_surf_alb_dir, i)     &
     &            +flandg_g(l) * (weight_690nm(i)                       &
     &            *land_albedo(list(l), 1)                              &
     &            +(1.0-weight_690nm(i))*land_albedo(list(l), 3))
              else
                rho_alb(l, ip_surf_alb_diff, i)                         &
     &            =weight_690nm(i)*land_albedo(list(l), 2)              &
     &            +(1.0-weight_690nm(i))*land_albedo(list(l), 4)
                rho_alb(l, ip_surf_alb_dir, i)                          &
     &            =weight_690nm(i)*land_albedo(list(l), 1)              &
     &            +(1.0-weight_690nm(i))*land_albedo(list(l), 3)
              endif
!
            else
!
!             If MOSES-II is not on, we cannot have coastal tiling,
!             so this point must be completely land.
              rho_alb(l, ip_surf_alb_diff, i)                           &
     &          =land_alb(list(l))
              rho_alb(l, ip_surf_alb_dir, i)                            &
     &          =land_alb(list(l))
!
            endif
!
          endif
!
        enddo


      enddo
!
!
!
!     Set the surface basis functions for a Lambertian surface.
      n_brdf_basis_fnc=1
!     By setting F_{1,0,0,0} equal to 4 we can set rho_alb equal to
!     the diffuse albedo.
      f_brdf(1, 0, 0, 0)=4.0
      if (ls_brdf_trunc /= 0) then
        write(iu_err, '(/a)')                                           &
     &    'error: the order of surface truncation is too high.'
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
        if (l_MOSES_II .AND. l_ctile) then
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
          do ll=1, nlit
            l=list(ll)
            if ( (flandg(l) < 1.0) .AND.                                &
     &           ( (flandg(l) > 0.0) .OR.                               &
     &             ( (ice_fraction(l) > 0.0) .AND.                      &
     &               (ice_fraction(l) < 1.0) ) ) ) then
              n_point_tile=n_point_tile+1
              list_tile(n_point_tile)=ll
              list_tile_outer(n_point_tile)=l
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
          do ll=1, nlit
            l=list(ll)
            if ( (.not.land(l)).and.                                    &
     &           (ice_fraction(l) >  0.0).and.                          &
     &           (ice_fraction(l) <  1.0) ) then
              n_point_tile=n_point_tile+1
              list_tile(n_point_tile)=ll
              list_tile_outer(n_point_tile)=l
            endif
          enddo
!
        endif
!
!       Now assign the tiled surface properties at points where tiling
!       is active.
        do i=1, n_band
!
!         The oceanic surface.
          rho_alb_tile(1:n_point_tile                                   &
     &      , ip_surf_alb_dir, ip_ocean_tile, i)                        &
     &      =open_sea_albedo(list_tile_outer(1:n_point_tile), 1)
          rho_alb_tile(1:n_point_tile                                   &
     &      , ip_surf_alb_diff, ip_ocean_tile, i)                       &
     &      =open_sea_albedo(list_tile_outer(1:n_point_tile), 2)
!
          if (l_MOSES_II) then
!
            if (l_ctile) then
!
!             With coastal tiling there is a real distinction between
!             land and seaice.
!
!             Seaice
              rho_alb_tile(1:n_point_tile                               &
     &          , ip_surf_alb_dir, ip_seaice_tile, i)                   &
     &          =sice_alb(list_tile_outer(1:n_point_tile))
              rho_alb_tile(1:n_point_tile                               &
     &          , ip_surf_alb_diff, ip_seaice_tile, i)                  &
     &          =sice_alb(list_tile_outer(1:n_point_tile))
!
!             Land
              rho_alb_tile(1:n_point_tile                               &
     &          , ip_surf_alb_dir, ip_land_tile, i)                     &
     &          =weight_690nm(i)                                        &
     &          *land_albedo(list_tile_outer(1:n_point_tile), 1)        &
     &          +(1.0-weight_690nm(i))                                  &
     &          *land_albedo(list_tile_outer(1:n_point_tile), 3)
              rho_alb_tile(1:n_point_tile                               &
     &          , ip_surf_alb_diff, ip_land_tile, i)                    &
     &          =weight_690nm(i)                                        &
     &          *land_albedo(list_tile_outer(1:n_point_tile), 2)        &
     &          +(1.0-weight_690nm(i))                                  &
     &          *land_albedo(list_tile_outer(1:n_point_tile), 4)
            else
!
!             The land albedo fields contain the valuse for sea-ice.
              rho_alb_tile(1:n_point_tile                               &
     &          , ip_surf_alb_dir, ip_seaice_tile, i)                   &
     &          =weight_690nm(i)                                        &
     &          *land_albedo(list_tile_outer(1:n_point_tile), 1)        &
     &          +(1.0-weight_690nm(i))                                  &
     &          *land_albedo(list_tile_outer(1:n_point_tile), 3)
              rho_alb_tile(1:n_point_tile                               &
     &          , ip_surf_alb_diff, ip_seaice_tile, i)                  &
     &          =weight_690nm(i)                                        &
     &          *land_albedo(list_tile_outer(1:n_point_tile), 2)        &
     &          +(1.0-weight_690nm(i))                                  &
     &          *land_albedo(list_tile_outer(1:n_point_tile), 4)
!
            endif
!
          else
!
!           The land albedo field contains the data for seaice.
            rho_alb_tile(1:n_point_tile                                 &
     &        , ip_surf_alb_dir, ip_seaice_tile, i)                     &
     &        =land_alb(list_tile_outer(1:n_point_tile))
            rho_alb_tile(1:n_point_tile                                 &
     &        , ip_surf_alb_diff, ip_seaice_tile, i)                    &
     &        =land_alb(list_tile_outer(1:n_point_tile))
!
          endif
!
        enddo
      endif
!
!
!
      return
      END SUBROUTINE r2_set_surface_field_sw
#endif
