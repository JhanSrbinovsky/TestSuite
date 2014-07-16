#if defined(A70_1Z)
#if defined(A01_3Z) ||  defined(A02_3Z)
! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************
!
!+ subroutine to set the mixing ratios of gases.
!
! Purpose:
!   the full array of mass mixing ratios of gases is filled.
!
! Method:
!   the arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. for well-mixed
!   gases the constant mixing ratios are fed into this array.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
!
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- --------------------------------------------------------------------
!+ subroutine to set thermodynamic properties
!
! Purpose:
!   pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.
!
! Method:
!   straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to assign properties of clouds.
!
! Purpose:
!   the fractions of different types of clouds and their microphysical
!   preoperties are set.
!
! Method:
!   straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
!
!
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   the parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
!
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to set fields of aerosols.
!
! Purpose:
!   the mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to set fields of climatological aerosols in hadcm3.
!
! Purpose:
!   this routine sets the mixing ratios of climatological aerosols.
!   a separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   the climatoogy used here is the one devised for hadcm3.
!
! Method:
!   straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
!
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      subroutine r2_set_aero_clim_hadcm3(n_profile, nlevs, n_layer      &
     &   , i_gather, l_extra_top                                        &
     &   , l_clim_aero_hgt, bl_depth, t, n_levels_bl                    &
     &   , land, lying_snow, pstar, p_layer_boundaries, trindx          &
     &   , aerosol_mix_ratio_clim                                       &
     &   , nd_field, nd_profile, nd_layer                               &
     &   )
!
!
!
      implicit none
!
!
!     comdecks included.
#include "c_g.h"
#include "c_r_cp.h"
!
!     dummy arguments.
!
!     sizes of arrays:
      integer                                                           &
                !, intent(in)
     &     nd_field                                                     &
!             field size in calling program
     &   , nd_profile                                                   &
!             size of array of profiles
     &   , nd_layer
!             maximum number of layers
!
!     actual sizes used:
      integer                                                           &
                !, intent(in)
     &     n_profile                                                    &
!             number of profiles
     &   , nlevs                                                        &
!             number of atmospheric layers used outside the radiation
!             scheme
     &   , n_layer
!             number of atmospheric layers seen in radiation
!
!     variables related to options for setting the field
      logical, intent(in) ::                                            &
     &     l_clim_aero_hgt                                              &
!             flag to use the depth of the boundary layer to set
!             the layers in which a boundary-layer aerosol is used
     &   , l_extra_top
!             flag to use an extra top layer in radiative calculations
      integer, intent(in) ::                                            &
     &     n_levels_bl
!             common number of layers taken to be occupied by the
!             boundary-layer aerosol if this the boundary layer
!             depth is not used to determine the number separately
!             at each grid-point
      real, intent(in) ::                                               &
     &     bl_depth(nd_field)                                           &
!             depth of the boundary layer
     &   , t(nd_profile, nd_layer)
!             temperatures of atmospheric layers
!
!     gathering array:
      integer                                                           &
                !, intent(in)
     &     i_gather(nd_field)
!             list of points to gather
!
!     general atmospheric properties:
      integer                                                           &
                !, intent(in)
     &     trindx(nd_field)
!             layer boundary of tropopause
      real                                                              &
                !, intent(in)
     &     pstar(nd_field)                                              &
!             surface pressures
     &,    p_layer_boundaries(nd_field,0:nlevs)
!             pressure at boundaries of layers
!
!     surface fields
      logical                                                           &
                !, intent(in)
     &     land(nd_field)
!             land-sea mask
      real                                                              &
                !, intent(in)
     &     lying_snow(nd_field)
!             depth of lying snow
!
      real                                                              &
                !, intent(out)
     &     aerosol_mix_ratio_clim(nd_profile, nd_layer, 5)
!             mixing ratios of climatological aerosols
!
!
!
!     local variables:
      integer                                                           &
     &     i                                                            &
!             loop variable
     &   , i_um                                                         &
!             index of a level in the um's upward convention
     &   , j                                                            &
!             loop variable
     &   , l                                                            &
!             loop variable
     &   , lg                                                           &
!             index for gathering
     &   , bl_top_lyr(nd_profile)                                       &
!             topmost layer occupied by boundary layer aerosol,
!             indexed using the top-down convention of radiation.
     &   , n_points_set_z
!             number of points where the top of the boundary layer
!             has not yet been reached
      logical                                                           &
     &     l_set_z(nd_profile)
!             array to flag points where the top of the boundary layer
!             has not yet been reached
      real                                                              &
     &     pressure_wt(nd_field)                                        &
!             array for scaling aerosol amounts for different surface
!             pressures
     &   , p_toa(nd_profile)                                            &
!             pressure at the top of the atmosphere seen by radiation
     &   , z_bottom(nd_profile)                                         &
!             height of the bottom of the current layer above
!             the surface
     &   , dz
!             depth of the current layer
!
!     total column mass (kg m-2) of each aerosol species in
!     the boundary layer, the free troposphere and the stratosphere
!     respectively. this model assumes that there are five aerosols.
      real                                                              &
     &     bl_oceanmass(5)                                              &
     &   , bl_landmass(5)                                               &
     &   , freetrop_mass(5)                                             &
     &   , strat_mass(5)
!
!     initialization for the climatological aerosol model
      data bl_landmass/2.77579e-5, 6.70018e-5, 0.0, 9.57169e-7, 0.0/
      data bl_oceanmass/1.07535e-5, 0.0, 2.043167e-4, 0.0, 0.0/
      data freetrop_mass/3.46974e-6, 8.37523e-6, 0.0, 1.19646e-7, 0.0/
      data strat_mass/0.0, 0.0, 0.0, 0.0, 1.86604e-6/
!
!
!
!     tropospheric aerosol loading is a simple function of surface
!     pressure: halving pstar halves the tropospheric aerosol burden.
!     the stratospheric burden is independent of pstar.  note the
!     factor multipling aerosol amounts uses a reference pressure
!     of 1013 mbars.
      do l=1, n_profile
        pressure_wt(l)=pstar(i_gather(l))*(1.0/1.013e5)
      end do
!
!     for each of the 5 aerosol species, the column amount in the
!     boundary layer, free troposphere and stratosphere are known for
!     a standard atmosphere over ocean and land. these can be used
!     to find mixing ratios for the um by dividing total aerosol by
!     total air mass (and using pressure weighting in the
!     troposphere).
!
!     firstly, mixing ratios are set for the 5 aerosol species in the
!     stratosphere.
      if (l_extra_top) then
!       with an extra layer the radiative atmosphere is notionally
!       extended to zero pressure.
        do l=1, n_profile
          p_toa(l)=0.0e+00
        enddo
      else
!       otherwise the top of the atmosphere seen elsewhere in the
!       model is used.
        do l=1, n_profile
          p_toa(l)=p_layer_boundaries(i_gather(l), nlevs)
        enddo
      endif
      do i=1,5
        do l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio_clim(l, n_layer+1-trindx(lg), i)            &
     &      =strat_mass(i)*g/                                           &
     &        (p_layer_boundaries(lg,trindx(lg)-1)                      &
     &       - p_toa(l))
        end do
      end do
      do i=1,5
        do l=1, n_profile
          lg=i_gather(l)
            do j=(trindx(lg)+1),n_layer
              aerosol_mix_ratio_clim(l, n_layer+1-j, i)=                &
     &          aerosol_mix_ratio_clim(l, n_layer+1-trindx(lg), i)
            end do
         end do
       end do
!
!      at each point the number of layers considered to contain
!      boundary layer aerosol must be determined. in the original
!      form of the scheme this number is fixed, but in the newer
!      form it is determined from the boundary layer depth.
       if (l_clim_aero_hgt) then
!
!        initialize:
         do l=1, n_profile
           bl_top_lyr(l)=n_layer
           l_set_z(l)=.true.
           z_bottom(l)=0.0e+00
         enddo
         n_points_set_z=n_profile
           i=n_layer
!
         do while (n_points_set_z >  0)
!
!          assign the um's indexing over layers: the um now indexes the
!          boundaries of layers from 0 at the bottom to nlevs at the
!          top, while the radiation code uses the same range but starts
!          at the top (possibly with an extra layer at the top of
!          the model). i and i_um refer to layers, not to the edges of
!          layers. bl_top_lyr now holds the topmost layer containing
!          boundary layer aerosol indexed as in the radiation code.
           i_um=n_layer+1-i
           do l=1, n_profile
             if (l_set_z(l)) then
               lg=i_gather(l)
               dz=r*t(l, i)*log(p_layer_boundaries(lg, i_um-1)          &
     &           /p_layer_boundaries(lg, i_um))/g
               if ( (max(bl_depth(lg), 5.0e+02) >=                      &
     &               z_bottom(l)+0.5e+00*dz).and.                       &
     &              (i_um <= trindx(lg)-2) ) then
!                the top of the boundary layer is above the middle
!                of this layer, so we take the layer to contain
!                boundary layer aerosol, provisionally setting
!                the index of the top of the boundary layer to this
!                level: it will be overwritten if higher layers are
!                also in the boundary layer. a upper limit is applied
!                to the number of the layers that can be filled with
!                the boundary layer aerosol so that there is at least
!                one layer in the free troposphere.
                 bl_top_lyr(l)=i
!                increment the height of the bottom of the layer
!                for the next pass.
                 z_bottom(l)=z_bottom(l)+dz
               else
!                the top of the boundary layer has been passed at this
!                point: it does not need to be considered further.
                 l_set_z(l)=.false.
                 n_points_set_z=n_points_set_z-1
               endif
             endif
           enddo
           i=i-1
         enddo
       else
         do l=1, n_profile
!          we do not allow the stratosphere to meet the boundary
!          layer, so we ensure that, counting up from the surface,
!          the highest layer allowed to be filled with boundary-
!          layer aerosol is the trindx-2nd: as a result, there must
!          be at least one layer in the free troposphere.
           bl_top_lyr(l)=n_layer+1                                      &
     &       -min(n_levels_bl, trindx(i_gather(l))-2)
         enddo
       endif
!
!      now, the mixing ratios are set for the 5 aerosol species
!      in the free troposphere. initially set the mixing ratio
!      in the lowest layer of the free troposphere.
       do i=1,5
         do l=1, n_profile
           lg=i_gather(l)
           aerosol_mix_ratio_clim(l, bl_top_lyr(l)-1, i)                &
     &       =freetrop_mass(i)*g*                                       &
     &       pressure_wt(l)/                                            &
     &       (p_layer_boundaries(lg, n_layer+1-bl_top_lyr(l))-          &
     &        p_layer_boundaries(lg,trindx(lg)-1) )
         end do
       end do
!      fill in the remaining levels where necessary.
       do l=1, n_profile
         lg=i_gather(l)
         do i=1,5
           do j=n_layer+2-trindx(lg), bl_top_lyr(l)-2
             aerosol_mix_ratio_clim(l, j, i)=                           &
     &         aerosol_mix_ratio_clim(l, bl_top_lyr(l)-1, i)
           end do
         end do
       end do
!
!      now, the boundary layer mixing ratios are set for the
!      5 aerosol species. a continental aerosol is used over most land
!      areas, but not over ice sheets, which are identified by the
!      criterion used in the boundary layer scheme that the mass of
!      lying snow exceeds 5000 kgm-2. over ice sheets a maritime
!      aerosol is used.
       do i=1,5
         do l=1, n_profile
           lg=i_gather(l)
           if ( land(lg).and.(lying_snow(lg) <  5.0e+03) ) then
             aerosol_mix_ratio_clim(l, bl_top_lyr(l), i)                &
     &         =bl_landmass(i)*g*pressure_wt(l)                         &
     &           /(pstar(lg)                                            &
     &           -p_layer_boundaries(lg, n_layer+1-bl_top_lyr(l)))
           else
             aerosol_mix_ratio_clim(l, bl_top_lyr(l), i)                &
     &         =bl_oceanmass(i)*g*pressure_wt(l)                        &
     &           /(pstar(lg)                                            &
     &           -p_layer_boundaries(lg, n_layer+1-bl_top_lyr(l)))
           end if
         end do
       end do
       do i=1,5
         do l=1, n_profile
           do j=bl_top_lyr(l)+1, n_layer
            aerosol_mix_ratio_clim(l,j,i)=                              &
     &         aerosol_mix_ratio_clim(l, bl_top_lyr(l), i)
           end do
         end do
       end do
!
!
!
      return
      END SUBROUTINE r2_set_aero_clim_hadcm3
!+ subroutine to calculate the total cloud cover.
!
! Purpose:
!   the total cloud cover at all grid-points is determined.
!
! Method:
!   a separate calculation is made for each different assumption about
!   the overlap.
!
! Current owner of code: J.-C. Thelen
!
! History:
!
!     Version            Date                   Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to implement the mrf umist parametrization.
!
! Purpose:
!   effective radii are calculated in accordance with this
!   parametrization.
!
! Method:
!   the number density of ccn is found from the concentration
!   of aerosols, if available. this yields the number density of
!   droplets: if aerosols are not present, the number of droplets
!   is fixed. effective radii are calculated from the number of
!   droplets and the lwc. limits are applied to these values. in
!   deep convective clouds fixed values are assumed.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to set the actual process options for the radiation code.
!
! Purpose:
!   to set a consistent set of process options for the radiation.
!
! Method:
!   the global options for the spectral region are compared with the
!   contents of the spectral file. the global options should be set
!   to reflect the capabilities of the code enabled in the model.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
#endif
#endif
