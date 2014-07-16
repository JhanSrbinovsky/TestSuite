#if defined(A70_1Z)
#if defined(A01_3Z) ||  defined(A02_3Z)
! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************
!
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
      subroutine r2_set_aerosol_field(ierr                              &
     &   , n_profile, nlevs, n_layer, n_aerosol, type_aerosol           &
     &   , i_gather, l_extra_top                                        &
     &   , l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero       &
     &   , bl_depth, t, n_levels_bl, l_murk_rad, aero_meso              &
     &   , l_use_sulpc_direct                                           &
     &   , sulp_dim1, sulp_dim2                                         &
     &   , accum_sulphate, aitken_sulphate                              &
     &   , l_use_seasalt_direct, salt_dim_a, salt_dim_b                 &
     &   , sea_salt_film, sea_salt_jet, p                               &
     &   , l_use_soot_direct, soot_dim1, soot_dim2                      &
     &   , fresh_soot, aged_soot                                        &
     &   , l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic       &
     &   , n_arcl_species, n_arcl_compnts, i_arcl_compnts               &
     &   , l_use_arcl, arcl_dim1, arcl_dim2, arcl                       &
     &   , land, lying_snow, pstar                                      &
     &   , p_layer_boundaries                                           &
     &   , trindx                                                       &
     &   , aerosol_mix_ratio                                            &
     &   , nd_field, nd_profile, nd_layer, nd_aerosol_species           &
     &   )
!
!
!
      implicit none
!
!
!     comdecks included.
#include "c_g.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
#include "aercmp3a.h"
#include "arcl_ids.h"
#include "c_r_cp.h"
#include "c_pi.h"
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
     &   , nd_layer                                                     &
!             maximum number of layers
     &   , nd_aerosol_species
!             maximum number of aerosol species
!
      integer                                                           &
                !, intent(out)
     &     ierr
!             error flag
!
!     actual sizes used:
      integer                                                           &
                !, intent(in)
     &     n_profile                                                    &
!             number of profiles
     &   , nlevs                                                        &
!             number of layers used outside the radiation scheme
     &   , n_layer                                                      &
!             number of layers seen by radiation
     &   , n_levels_bl                                                  &
!             number of layers occupied by boundary-layer aerosol
!             if l_clim_aero_hgt is false.
     &   , n_aerosol                                                    &
!             number of aerosols in spectral file
     &   , type_aerosol(nd_aerosol_species)
!             actual types of aerosols
!
!     gathering array:
      integer                                                           &
                !, intent(in)
     &     i_gather(nd_field)
!             list of points to gather
!
!     flag for the climatological aerosol distribution.
      logical                                                           &
                   !, intent(in)
     &     l_climat_aerosol                                             &
!             flag for climatological aerosol distribution
     &   , l_clim_aero_hgt                                              &
!             flag to use the boundary layer depth in setting the
!             climatological aerosol
     &   , L_HadGEM1_Clim_Aero                                          &
!             Flag to use HadGEM1 setting for climatological aerosols
     &   , l_murk_rad                                                   &
!             flag for mesoscale model aerosol
     &   , l_extra_top
!             flag to include an extra top layer in the radiation
!             scheme
!
!     variables for the sulphur cycle:
      logical                                                           &
                   !, intent(in)
     &     l_use_sulpc_direct
!             flag to use sulphur cycle for direct effect
      integer                                                           &
                   !, intent(in)
     &     sulp_dim1,sulp_dim2
!             dimensions for _sulphate arrays, (p_field,p_levels or 1,1)
      real                                                              &
                !, intent(in)
     &     accum_sulphate(sulp_dim1, sulp_dim2)                         &
!             mass mixing ratios of accumulation mode aerosol
     &   , aitken_sulphate(sulp_dim1, sulp_dim2)
!             mass mixing ratios of aitken mode aerosol
!
! declare soot variables:
      logical l_use_soot_direct !use direct rad. effect of soot aerosol
      integer soot_dim1,soot_dim2
                !dimensions for soot arrays, (p_field,p_levels or 1,1)
      real fresh_soot(soot_dim1, soot_dim2)                             &
                                                 ! mmr of fresh soot
     &   , aged_soot(soot_dim1, soot_dim2)       ! mmr of aged soot
!
      logical                                                           &
     &     l_use_seasalt_direct
!             flag for direct effect of interactive sea-salt aerosol
!
      integer                                                           &
     &     salt_dim_a, salt_dim_b
!             array sizes of sea-salt aerosols
!                (either p_field,p_levels or 1,1)
!
      real                                                              &
     &     sea_salt_film(salt_dim_a, salt_dim_b)                        &
!             on input, number concentration (m-3) of film-mode
!             sea-salt aerosols; converted to mass mixing ratio.
     &   , sea_salt_jet(salt_dim_a, salt_dim_b)
!             on input, number concentration (m-3) of jet-mode
!             sea-salt aerosols; converted to mass mixing ratio.

! declare biogenic variables
      logical l_use_biogenic ! flag for biogenic aerosol direct effect
      ! dimensions of biogenic array
      integer biogenic_dim1, biogenic_dim2
      ! biogenic aerosol mass-mixing ratios
      real biogenic(biogenic_dim1, biogenic_dim2)

!
! Aerosol climatology for NWP:
!
#include "arcl_dim.h"

      ! Number of requested species within the climatology
      integer n_arcl_species
      
      ! Corresponding number of requested components
      integer n_arcl_compnts
      
      ! Model switches for each species
      logical, dimension(NPD_ARCL_SPECIES) :: l_use_arcl
      
      ! Array index of components
      integer, dimension(NPD_ARCL_COMPNTS) :: i_arcl_compnts
      
      ! Array dimensions
      integer arcl_dim1, arcl_dim2
      
      ! Mass-mixing ratios
      real                                                              &
     &    arcl(arcl_dim1, arcl_dim2, n_arcl_compnts)

      real                                                              &
     &     aero_meso(nd_field, nlevs)
!             mixing ratio of 'urban' aerosol of mesoscale model
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
!             land sea mask
      real                                                              &
                !, intent(in)
     &     lying_snow(nd_field)                                         &
!             depth of lying snow
     &   , bl_depth(nd_field)                                           &
!             depth of the boundary layer
     &   , t(nd_profile, nd_layer)                                      &
!             temperatures of atmospheric layers
     &   , p(nd_profile, nd_layer)
!             pressures at layer centres
!
      real                                                              &
                !, intent(out)
     &     aerosol_mix_ratio(nd_profile, nd_layer                       &
     &        , nd_aerosol_species)
!             mixing ratios of aerosols
!
!
!
!     local variables:
      integer                                                           &
     &     i                                                            &
!             loop variable
     &   , j                                                            &
!             loop variable
     &   , l                                                            &
!             loop variable
     &   , lg                                                           &
!             index for gathering
     &   , bltop                                                        &
!             index of upper boundary of planetary boundary layer
     &   , i_aerosol                                                    &
!             actual type of aerosol being considered
     &   , i_top_copy
!             topmost radiative layer to be set directly from the
!             profile supplied
!
!
!     arrays for the climatological aerosol model
      integer :: ii
!       variable used in initialization of arrays.
      logical, dimension(npd_aerosol_component) :: l_in_climat =        &
     &  (/ (.true., ii=1, 4), .false., .true.,                          &
     &     (.false., ii=7, npd_aerosol_component) /)
!       flags to indicate which aerosols are included in the
!       climatology: this may be used to enable various components
!       to be replaced by fully prognostic schemes.
      integer, dimension(npd_aerosol_component) :: i_clim_pointer =     &
     &  (/ (ii, ii=1, 4), 0, 5,                                         &
     &     (0, ii=7, npd_aerosol_component) /)
!       pointers to the indices of the original climatological
!       aerosol model.
      real                                                              &
     &     aerosol_mix_ratio_clim(nd_profile, nd_layer, 5)
!             mixing ratios of the climatological aerosols
!
      real                                                              &
     &     mode_radius_ss_film                                          &
!            mode radius of film-mode sea-salt aerosol (m)
     &   , mode_radius_ss_jet                                           &
!            mode radius of jet-mode sea-salt aerosol (m)
     &   , sigma_ss_film                                                &
!            geometric standard deviation of film-mode sea-salt aerosol
     &   , sigma_ss_jet                                                 &
!            geometric standard deviation of jet-mode sea-salt aerosol
     &   , density_sea_salt
!            bulk density of sodium chloride, taken as representative
!            of sea-salt aerosol (kg m-3)
!
      integer im
      real, dimension(npd_aerosol_component) :: meso_frac =             &
     &  (/ 0.61, 0.17, 0.0, 0.22, (0.0, im=5,npd_aerosol_component) /)
      real inv17

      parameter(                                                        &
     &     mode_radius_ss_film=0.1e-06                                  &
     &   , mode_radius_ss_jet=1.0e-06                                   &
     &   , sigma_ss_film=1.9                                            &
     &   , sigma_ss_jet=2.0                                             &
     &   , density_sea_salt=2160.0                                      &
     &   )
!
!
!
!     subroutines called:
      external                                                          &
     &     r2_set_aero_clim_hadcm3
!
!
!
!     if an extra layer is used in the radiation scheme, any
!     aerosol profiles supplied will be copied into the profiles
!     seen by the radiation scheme starting with the second layer
!     down, rather than the first.
      i_top_copy=1
      if (l_extra_top) i_top_copy=2
!
!  Use HadGEM1 settings for climatological aerosols.
!
      if (L_HadGEM1_Clim_Aero) then
        do ii=1,4
          L_IN_CLIMAT(ii) = .false.
        enddo
      endif
!
!
!  use climatological soot if climatological aerosols are on and not
!  using interactive soot.
       l_in_climat(ip_soot) = l_in_climat(ip_soot)                      &
     &                     .and.(.not.l_use_soot_direct)
!
!  use climatological "oceanic" aerosol if climatological aerosols have
!  been selected and not using sea-salt parametrization.
      l_in_climat(ip_oceanic) = l_in_climat(ip_oceanic)                 &
     &                          .and.(.not. l_use_seasalt_direct)
!
!  the climatological water-soluble aerosol should not be used
!  if the direct effects of sulphate aerosols are included.
!  (note that this differs the situation applying in earlier
!  versions of the model, such as hadam3).
      l_in_climat(ip_water_soluble) = l_in_climat(ip_water_soluble)     &
     &                     .and.(.not.l_use_sulpc_direct)
!
!
!
      if (l_climat_aerosol) then
!
!        set the mixing ratios of the climatological aerosols
!        used in the climatology of hadcm3. a separate subroutine
!        is used to ensure bit-reproducible results by using
!        earlier code. this could be altered if a new climatology were
!        used.
!
! DEPENDS ON: r2_set_aero_clim_hadcm3
         call r2_set_aero_clim_hadcm3(n_profile, nlevs, n_layer         &
     &      , i_gather, l_extra_top                                     &
     &      , l_clim_aero_hgt, bl_depth, t, n_levels_bl                 &
     &      , land, lying_snow, pstar, p_layer_boundaries, trindx       &
     &      , aerosol_mix_ratio_clim                                    &
     &      , nd_field, nd_profile, nd_layer                            &
     &      )
!
      endif
!
!
!     the aerosols required by for the calculation should have been
!     selected when the spectral file was read in. each type should
!     be set appropriately.
!
      do j=1, n_aerosol
!
         i_aerosol=type_aerosol(j)
!
         if (l_climat_aerosol.and.l_in_climat(i_aerosol)) then
!
!           here mxing ratios in all layers are set because the
!           possible presence of an extra layer has been allowed
!           for in the subroutine.
            do i=1, n_layer
               do l=1, n_profile
                  aerosol_mix_ratio(l, i, j)                            &
     &               =aerosol_mix_ratio_clim(l, i                       &
     &               , i_clim_pointer(i_aerosol))
               enddo
            enddo
!
! For the other aerosol species, they may be supplied from dedicated
! schemes (e.g. climate model) or climatologies (e.g. NWP model).
! The former are indicated by L_USE_<SPECIES>, the latter by
! L_USE_ARCL(IP_ARCL_<SPECIES>). Should both logical be set to true
! for a given aerosol species (i.e. an aerosol is provided by a 
! dedicated scheme AND a climatology), the climatology wins and will 
! be seen by radiation.
!
! Compared to version 3C, this version does not support all the
! interactive aerosol components. All aerosols from NWP climatology
! are supported, however.
!
!
         else if (i_aerosol == ip_accum_sulphate) then
!
!           aerosols related to the sulphur cycle (note that dissolved
!           sulphate does not contribute to the direct effect):
!
            if (l_use_arcl(IP_ARCL_SULP)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)                          &
                       =arcl(lg, n_layer+1-i,                           &
     &                       i_arcl_compnts(IP_ARCL_SULP_AC))
                 end do
               end do
             else if (l_use_sulpc_direct) then
               do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)                          &
     &                 =accum_sulphate(lg, n_layer+1-i)
                 end do
              end do
            end if
!
         else if (i_aerosol == ip_aitken_sulphate) then
            
            if (l_use_arcl(IP_ARCL_SULP)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)                          &
                       =arcl(lg, n_layer+1-i,                           &
     &                       i_arcl_compnts(IP_ARCL_SULP_AK))
                 end do
               end do
             else if (l_use_sulpc_direct) then
               do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)                          &
     &                 =aitken_sulphate(lg, n_layer+1-i)
                 end do
              end do
            end if
!
         else if (i_aerosol == ip_fresh_soot) then
            
            if (l_use_arcl(IP_ARCL_BLCK)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)                          &
                       =arcl(lg, n_layer+1-i,                           &
     &                       i_arcl_compnts(IP_ARCL_BLCK_FR))
                 end do
               end do
             else if (l_use_soot_direct) then
               do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)                          &
     &                 =fresh_soot(lg, n_layer+1-i)
                 end do
              end do
            end if
!
         else if (i_aerosol == ip_aged_soot) then
            
            if (l_use_arcl(IP_ARCL_BLCK)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)                          &
                       =arcl(lg, n_layer+1-i,                           &
     &                       i_arcl_compnts(IP_ARCL_BLCK_AG))
                 end do
               end do
             else if (l_use_soot_direct) then
               do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)                          &
     &                 =aged_soot(lg, n_layer+1-i)
                 end do
              end do
            end if
!
         else if (i_aerosol == ip_seasalt_film) then
            
            if (l_use_arcl(IP_ARCL_SSLT)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                          (arcl(lg, n_layer+1-i,                        &
     &                          i_arcl_compnts(IP_ARCL_SSLT_FI))*4.0*pi &
     &                   *density_sea_salt*(mode_radius_ss_film**3.0)   &
     &                   *exp(4.5*(alog(sigma_ss_film))**2.0))          &
     &                   /(3.0*p(l, i)/(r*t(l, i)))
                 end do
              end do
            else if (l_use_seasalt_direct) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
     &                    (sea_salt_film(lg, n_layer+1-i)*4.0*pi        &
     &                   *density_sea_salt*(mode_radius_ss_film**3.0)   &
     &                   *exp(4.5*(alog(sigma_ss_film))**2.0))          &
     &                   /(3.0*p(l, i)/(r*t(l, i)))
                 end do
              end do
            end if
!
         else if (i_aerosol == ip_seasalt_jet) then

            if (l_use_arcl(IP_ARCL_SSLT)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                          (arcl(lg, n_layer+1-i,                        &
     &                          i_arcl_compnts(IP_ARCL_SSLT_JT))*4.0*pi &
     &                   *density_sea_salt*(mode_radius_ss_jet**3.0)    &
     &                   *exp(4.5*(alog(sigma_ss_jet))**2.0))           &
     &                   /(3.0*p(l, i)/(r*t(l, i)))
                 end do
              end do
            else if (l_use_seasalt_direct) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
     &                    (sea_salt_jet(lg, n_layer+1-i)*4.0*pi         &
     &                   *density_sea_salt*(mode_radius_ss_jet**3.0)    &
     &                   *exp(4.5*(alog(sigma_ss_jet))**2.0))           &
     &                   /(3.0*p(l, i)/(r*t(l, i)))
                 end do
              end do
            end if
!
         else if (i_aerosol == ip_dust_1) then
         
            if (l_use_arcl(IP_ARCL_DUST)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_DUST_B1))
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_dust_2) then
         
            if (l_use_arcl(IP_ARCL_DUST)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_DUST_B2))
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_dust_3) then
         
            if (l_use_arcl(IP_ARCL_DUST)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_DUST_B3))
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_dust_4) then
         
            if (l_use_arcl(IP_ARCL_DUST)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_DUST_B4))
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_dust_5) then
         
            if (l_use_arcl(IP_ARCL_DUST)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_DUST_B5))
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_dust_6) then
         
            if (l_use_arcl(IP_ARCL_DUST)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_DUST_B6))
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_biomass_1) then
         
            if (l_use_arcl(IP_ARCL_BIOM)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_BIOM_FR))
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_biomass_2) then
         
            if (l_use_arcl(IP_ARCL_BIOM)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_BIOM_AG))
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_biogenic) then
         
            if (l_use_biogenic) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           biogenic(lg, n_layer+1-i)
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_ocff_fresh) then
         
            if (l_use_arcl(IP_ARCL_OCFF)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_OCFF_FR))
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_ocff_aged) then
         
            if (l_use_arcl(IP_ARCL_OCFF)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_OCFF_AG))
                 end do
              end do
              
            end if
!
         else if (i_aerosol == ip_delta) then
         
            if (l_use_arcl(IP_ARCL_DLTA)) then
              do i=i_top_copy, n_layer
                 do l=1, n_profile
                    lg=i_gather(l)
                    aerosol_mix_ratio(l, i, j)=                         &
                           arcl(lg, n_layer+1-i,                        &
     &                       i_arcl_compnts(IP_ARCL_DLTA_DL))
                 end do
              end do
              
            end if
!
         else
!
!           the options to the radiation code do not require this
!           aerosol to be considered: its mixing ratio is set to 0.
!           this block of code should not normally be executed,
!           but may be required for ease of including modifications.
!
            do i=1, n_layer
               do l=1, n_profile
                  lg=i_gather(l)
                  aerosol_mix_ratio(l, i, j)=0.0e+00
               enddo
            enddo
!
!
         endif
!
!        if using an extra top layer extrapolate the mixing ratio
!        from the adjacent layer, except in the case of climatological
!        aerosols which are specifically set.
         if (l_extra_top) then
           if ( ( (i_aerosol == ip_accum_sulphate).and.                 &
     &            (l_use_sulpc_direct.or.l_use_arcl(IP_ARCL_SULP)) )    &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_aitken_sulphate).and.                &
     &            (l_use_sulpc_direct.or.l_use_arcl(IP_ARCL_SULP)) )    &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_fresh_soot).and.                     &
     &            (l_use_soot_direct.or.l_use_arcl(IP_ARCL_BLCK)) )     &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_aged_soot).and.                      &
     &            (l_use_soot_direct.or.l_use_arcl(IP_ARCL_BLCK)) )     &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_seasalt_film).and.                   &
     &            (l_use_seasalt_direct.or.l_use_arcl(IP_ARCL_SSLT)))   &     
     &          .or.                                                    &
     &          ( (i_aerosol == ip_seasalt_jet).and.                    &
     &            (l_use_seasalt_direct.or.l_use_arcl(IP_ARCL_SSLT)))   &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_dust_1).and.                         &
     &             l_use_arcl(IP_ARCL_DUST))                            & 
     &          .or.                                                    &
     &          ( (i_aerosol == ip_dust_2).and.                         &
     &             l_use_arcl(IP_ARCL_DUST))                            &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_dust_3).and.                         &
     &             l_use_arcl(IP_ARCL_DUST))                            &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_dust_4).and.                         &
     &             l_use_arcl(IP_ARCL_DUST))                            &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_dust_5).and.                         &
     &             l_use_arcl(IP_ARCL_DUST))                            &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_dust_6).and.                         &
     &             l_use_arcl(IP_ARCL_DUST))                            &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_biomass_1).and.                      &
     &             l_use_arcl(IP_ARCL_BIOM))                            &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_biomass_2).and.                      &
     &             l_use_arcl(IP_ARCL_BIOM))                            &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_biogenic).and.                       &
     &             l_use_biogenic)                                      &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_ocff_fresh).and.                     &
     &             l_use_arcl(IP_ARCL_OCFF))                            &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_ocff_aged).and.                      &
     &             l_use_arcl(IP_ARCL_OCFF))                            &
     &          .or.                                                    &
     &          ( (i_aerosol == ip_delta).and.                          &
     &             l_use_arcl(IP_ARCL_DLTA))                            &
     &        ) then
             do l=1, n_profile
               aerosol_mix_ratio(l, 1, j)                               &
     &           =aerosol_mix_ratio(l, 2, j)
             enddo
           endif
         endif
!
      enddo
!
! only use aerosols in the boundary layer for the mesoscale model
      if (l_murk_rad) then
        inv17=1./1.7
        do j=1, n_aerosol
          do i=(nlevs-n_levels_bl+1),nlevs
            do l=1, n_profile
              lg=i_gather(l)
! note that aerosol mass mixing ratios are in units of 1e9 kg/kg
              aerosol_mix_ratio(l,i,j) = aerosol_mix_ratio(l,i,j)       &
     &          + aero_meso(lg,nlevs+1-i)*1e-9*meso_frac(j)*inv17
            enddo
          enddo
        enddo
      endif
!
!
      return
      END SUBROUTINE r2_set_aerosol_field
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
