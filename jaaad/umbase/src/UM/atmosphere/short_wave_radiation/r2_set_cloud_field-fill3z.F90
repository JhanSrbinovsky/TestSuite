#if defined(A70_1Z)
#if defined(A01_3Z) ||  defined(A02_3Z)
! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************
!
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
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      subroutine r2_set_cloud_field(n_profile, nlevs, n_layer, nclds    &
     &   , i_gather                                                     &
     &   , p, t, d_mass                                                 &
     &   , ccb, cct, cca, cccwp                                         &
     &   , lccwc1, lccwc2, lca_area, lca_bulk                           &
     &   , l_microphysics, l_aerosol_ccn                                &
     &   , sea_salt_film, sea_salt_jet                                  &
     &   , l_seasalt_ccn, salt_dim_a, salt_dim_b                        &
     &   , sulp_dim1, sulp_dim2, accum_sulphate, diss_sulphate          &
     &   , aitken_sulphate, lying_snow                                  &
     &   , l_cloud_water_partition, land_g, flandg_g                    &
     &   , i_cloud_representation, i_condensed_param                    &
     &   , condensed_min_dim, condensed_max_dim                         &
     &   , n_condensed, type_condensed                                  &
     &   , w_cloud, n_cloud_type, frac_cloud, l_local_cnv_partition     &
     &   , condensed_mix_rat_area, condensed_dim_char                   &
     &   , re_conv, re_conv_flag, re_strat, re_strat_flag               &
     &   , wgt_conv, wgt_conv_flag, wgt_strat, wgt_strat_flag           &
     &   , lwp_strat, lwp_strat_flag                                    &
     &   , ntot_diag, ntot_diag_flag                                    &
     &   , strat_lwc_diag, strat_lwc_diag_flag                          &
     &   , so4_ccn_diag, so4_ccn_diag_flag                              &
     &   , cond_samp_wgt, cond_samp_wgt_flag                            &
     &   , col_list, row_list, row_length, rows                         &
     &   , nd_field, nd_profile, nd_layer, id_ct, nd_aerosol_species    &
     &   , nd_cloud_component, nd_cloud_type                            &
     &   , n_cca_lev, Ntot_land, Ntot_sea                               &
     &   )
 
      Use cv_cntl_mod,                                                  &
          Only: lcv_3d_cca


 
      implicit none
!
!
!     comdecks included.
#include "cloud_component_pcf3z.h"
#include "cloud_type_pcf3z.h"
#include "cloud_representation_pcf3z.h"
#include "ice_cloud_param_pcf3z.h"
#include "c_0_dg_c.h"
#include "c_r_cp.h"
!
!
!     dimensions of arrays:
      integer, intent(in) :: row_length
!                              number of grid-points in ew-direction
!                              in the local domain
      integer, intent(in) :: rows
!                              number of grid-points in ns-direction
!                              in the local domain
      integer                                                           &
                !, intent(in)
     &     nd_field                                                     &
!             field size in calling program
     &   , nd_profile                                                   &
!             size of array of profiles
     &   , nd_layer                                                     &
!             maximum number of layers
     &   , nd_aerosol_species                                           &
!             maximum number of aerosol_species
     &   , nd_cloud_component                                           &
!             Size allocated for components in a cloud
     &   , nd_cloud_type                                                &
!             Size allocated for number of types of cloud
     &   , id_ct                                                        &
!             Topmost declared cloudy layer
     &   , sulp_dim1                                                    &
!             1st dimension of arrays of sulphate
     &   , sulp_dim2                                                    &
!             2nd dimension of arrays of sulphate
     &   , salt_dim_a                                                   &
!             1st dimension of arrays of sea-salt
     &   , salt_dim_b                                                   &
!             2nd dimension of arrays of sea-salt
     &   , n_cca_lev
!             number of levels for convective cloud amount
!
!     actual sizes used:
      integer                                                           &
                !, intent(in)
     &     n_profile                                                    &
!             number of profiles
     &   , nlevs                                                        &
!             number of layers used outside the radiation scheme
     &   , n_layer                                                      &
!             number of layers seen by the radiation code
     &   , nclds
!             number of cloudy levels
!
!     gathering array:
      integer                                                           &
                !, intent(in)
     &     i_gather(nd_field)
!             list of points to be gathered
      integer, intent(in) :: col_list(nd_field)
!                              ew indices of gathered points in the 2-d
!                              domain
      integer, intent(in) :: row_list(nd_field)
!                              ns indices of gathered points in the 2-d
!                              domain
!
!     thermodynamic fields:
      real                                                              &
                !, intent(in)
     &     p(nd_profile, nd_layer)                                      &
!             pressures
     &   , t(nd_profile, nd_layer)                                      &
!             temperatures
     &   , d_mass(nd_profile, nd_layer)
!             mass thicknesses of layers
!
!     convective clouds:
      integer                                                           &
                !, intent(in)
     &     ccb(nd_field)                                                &
!             base of convective cloud
     &   , cct(nd_field)
!             top of convective cloud
      real                                                              &
                !, intent(in)
     &     cca(nd_field,n_cca_lev)                                      &
!             fraction of convective cloud
     &   , cccwp(nd_field)
!             water path of convective cloud
      logical                                                           &
                !, intent(in)
     &     l_local_cnv_partition                                        &
!             flag to carry out the partitioning between ice
!             and water in convective clouds as a function of
!             the local temperature
     &   , l_seasalt_ccn
!              flag for sea-salt parametrization for ccn
!
!     layer clouds:
      real                                                              &
                !, intent(in)
     &     lccwc1(nd_field, nclds+1/(nclds+1))                          &
!             liquid water contents
     &   , lccwc2(nd_field, nclds+1/(nclds+1))                          &
!             ice water contents
     &   , lca_area(nd_field, nclds+1/(nclds+1))                        &
!             area coverage fractions of layer clouds
     &   , lca_bulk(nd_field, nclds+1/(nclds+1))
!             bulk coverage fractions of layer clouds
!
!     arrays for microphysics:
      logical                                                           &
                !, intent(in)
     &     l_microphysics                                               &
!             microphysical flag
     &   , l_aerosol_ccn                                                &
!             flag to use aerosols to find ccn
     &   , l_cloud_water_partition                                      &
!             flag to use prognostic cloud ice contents
     &   , land_g(nd_profile)
!             flag for land points
      real                                                              &
                !, intent(in)
     &     accum_sulphate(sulp_dim1, sulp_dim2)                         &
!             mixing ratios of accumulation-mode sulphate
     &   , aitken_sulphate(sulp_dim1, sulp_dim2)                        &
!             mixing ratios of aitken-mode sulphate
     &   , diss_sulphate(sulp_dim1, sulp_dim2)                          &
!             mixing ratios of dissolved sulphate
     &   , sea_salt_film(salt_dim_a, salt_dim_b)                        &
!             number concentration of film-mode sea-salt aerosol
     &   , sea_salt_jet(salt_dim_a, salt_dim_b)
!             number concentration of jet-mode sea-salt aerosol
!
!     representation of clouds
      integer                                                           &
                !, intent(in)
     &     i_cloud_representation
!             representation of clouds
!
!     parametrizations for clouds:
      integer                                                           &
                !, intent(in)
     &     i_condensed_param(nd_cloud_component)
!             types of parametrization used for condensed
!             components in clouds
!     limits on sizes of particles
      real                                                              &
                !, intent(in)
     &     condensed_min_dim(nd_cloud_component)                        &
!             minimum dimension of each condensed component
     &   , condensed_max_dim(nd_cloud_component)
!             maximum dimension of each condensed component
!
      real                                                              &
                !, Intent(IN)
     &     Ntot_land                                                    &
!             number of droplets over land / m-3
     &   , Ntot_sea
!             number of droplets over sea / m-3
!
!     assigned cloud fields:
      integer                                                           &
                !, intent(out)
     &     n_condensed                                                  &
!             number of condensed components
     &   , type_condensed(nd_cloud_component)                           &
!             types of condensed components
     &   , n_cloud_type
!             Number of types of clouds
      real                                                              &
                !, intent(out)
     &     w_cloud(nd_profile, id_ct: nd_layer)                         &
!             total amounts of cloud
     &   , frac_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)       &
!             fraction of each type of cloud
     &   , condensed_dim_char(nd_profile, id_ct: nd_layer               &
     &      , nd_cloud_component)                                       &
!             characteristic dimensions of cloudy components
     &   , condensed_mix_rat_area(nd_profile, id_ct: nd_layer           &
     &      , nd_cloud_component)                                       &
!             mass mixing ratios of condensed components using area cld
     &   , ntot_diag_g(nd_profile, nd_layer)                            &
!             diagnostic array for ntot (gathered)
     &   , strat_lwc_diag_g(nd_profile, nd_layer)                       &
!             diagnostic array for stratiform lwc (gathered)
     &   , so4_ccn_diag_g(nd_profile, nd_layer)
!             diagnostic array for so4 ccn mass conc (gathered)
!
!
      real                                                              &
     &     lying_snow(nd_field)                                         &
!            snow depth (>5000m = land ice sheet)
     &   , lying_snow_g(nd_profile)                                     &
!            gathered version of the above
     &   , flandg_g(nd_profile)
!            gathered global land fraction field
!
!     microphysical diagnostics:
      logical                                                           &
     &     re_conv_flag                                                 &
!             diagnose effective radius*weight for convective cloud
     &   , re_strat_flag                                                &
!             diagnose effective radius*weight for stratiform cloud
     &   , wgt_conv_flag                                                &
!             diagnose weight for convective cloud
     &   , wgt_strat_flag                                               &
!             diagnose weight for stratiform cloud
     &   , lwp_strat_flag                                               &
!             diagnose liquid water path*weight for stratiform cloud
     &   , ntot_diag_flag                                               &
!             diagnose droplet concentration*weight
     &   , strat_lwc_diag_flag                                          &
!             diagnose stratiform lwc*weight
     &   , so4_ccn_diag_flag                                            &
!             diagnose so4 ccn mass conc*cond. samp. weight
     &   , cond_samp_wgt_flag
!             diagnose conditional sampling weight
!
      real                                                              &
     &     re_conv(row_length, rows, nclds)                             &
!             effective radius*weight for convective cloud
     &   , re_strat(row_length, rows, nclds)                            &
!             effective radius*weight for stratiform cloud
     &   , wgt_conv(row_length, rows, nclds)                            &
!             weight for convective cloud
     &   , wgt_strat(row_length, rows, nclds)                           &
!             weight for stratiform cloud
     &   , lwp_strat(row_length, rows, nclds)                           &
!             liquid water path*weight for stratiform cloud
     &   , ntot_diag(row_length, rows, nclds)                           &
!             droplet concentration*weight
     &   , strat_lwc_diag(row_length, rows, nclds)                      &
!             stratiform lwc*weight
     &   , so4_ccn_diag(row_length, rows, nclds)                        &
!             so4 ccn mass conc*cond. samp. weight
     &   , cond_samp_wgt(row_length, rows, nclds)
!             conditional sampling weight
!
!
!
!     local variables:
      integer                                                           &
     &     i                                                            &
!             loop variable
     &   , l                                                            &
!             loop variable
     &   , lg
!             index to gather
!
      logical                                                           &
     &     l_glaciated_top(nd_profile)
!             logical for glaciated tops in convective cloud.
!
      real                                                              &
     &     liq_frac(nd_profile)                                         &
!             fraction of liquid cloud water
     &   , liq_frac_conv(nd_profile)                                    &
!             fraction of liquid water in convective cloud
     &   , t_gather(nd_profile)                                         &
!             gathered temperature for lsp_focwwil
     &   , total_mass(nd_profile)                                       &
!             total mass in convective cloud
     &   , cc_depth(nd_profile)                                         &
!             depth of convective cloud
     &   , condensed_mix_rat_bulk(nd_profile, id_ct: nd_layer           &
     &      , nd_cloud_component)                                       &
!             mass mixing ratios of condensed components using bulk cld
     &   , density_air(nd_profile, nd_layer)                            &
!             density of air
     &   , convective_cloud_layer(nd_profile)
!             amount of convective cloud in the current layer
!
!     Parameters for the aggregate ice scheme
      REAL, Parameter :: a0_agg_cold = 7.5094588E-04
      REAL, Parameter :: b0_agg_cold = 5.0830326E-07
      REAL, Parameter :: a0_agg_warm = 1.3505403E-04
      REAL, Parameter :: b0_agg_warm = 2.6517429E-05
      REAL, Parameter :: t_switch    = 216.208
      REAL, Parameter :: t0_agg      = 279.5
      REAL, Parameter :: s0_agg      = 0.05
!
!
!
!     set the components within the clouds. in the unified model we
!     have four components: stratiform ice and water and convective
!     ice and water.
      n_condensed=4
      type_condensed(1)=ip_clcmp_st_water
      type_condensed(2)=ip_clcmp_st_ice
      type_condensed(3)=ip_clcmp_cnv_water
      type_condensed(4)=ip_clcmp_cnv_ice
!
!
!
!     set the total amounts of cloud and the fractions comprised by
!     convective and stratiform components.
!
      if (i_cloud_representation == ip_cloud_conv_strat) then
!
         n_cloud_type=2
!
!  this cloud representation not available with new cloud microphysics
!
!        the clouds are divided into mixed-phase stratiform and
!        convective clouds: lsp_focwwil gives the partitioning between
!        ice and water in stratiform clouds and in convective cloud,
!        unless the option to partition as a function of the local
!        temperature is selected. within convective cloud the liquid
!        water content is distributed uniformly throughout the cloud.
!
!        convective cloud:
!
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)=0.0e+00
               condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)=0.0e+00
               condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)=0.0e+00
               condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)=0.0e+00
            enddo
         enddo
!
!
         if (l_local_cnv_partition) then
!
!           partition between ice and water using the relationships
!           given in bower et al. (1996, q.j. 122 p 1815-1844). ice
!           is allowed in a layer warmer than the freezing point
!           only if the top of the cloud is glaciated.
!
            do l=1, n_profile
!              min is required since cct may be 0 if there is no
!              convective cloud.
               l_glaciated_top(l)                                       &
     &            =(t(l, min(n_layer+2-cct(i_gather(l))                 &
     &            , n_layer-nclds+1)) <  tm)
            enddo

         else
!
!           partition between ice and water as directed by the
!           temperature in the middle of the top layer of the cloud.
!           the partitioning may be precalculated in this case.
!
            do l=1, n_profile
               t_gather(l)=t(l, min(n_layer+2-cct(i_gather(l))          &
     &            , n_layer-nclds+1))
            enddo
! DEPENDS ON: lsp_focwwil
            call lsp_focwwil(t_gather, n_profile, liq_frac_conv)
!
         endif
!
!
         do l=1, n_profile
            total_mass(l)=0.0e+00
         enddo
!
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               lg=i_gather(l)
               if ( (cct(lg) >= n_layer+2-i).and.                       &
     &              (ccb(lg) <= n_layer+1-i) ) then
                  total_mass(l)=total_mass(l)+d_mass(l, i)
               endif
            enddo
         enddo
!
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               lg=i_gather(l)
               if ( (cct(lg) >= n_layer+2-i).and.                       &
     &              (ccb(lg) <= n_layer+1-i) ) then
                  if (l_local_cnv_partition) then
!                    the partitioning is recalculated for each layer
!                    otherwise a generic value is used.
                     liq_frac_conv(l)=max(0.0e+00, min(1.0e+00          &
     &                  , 1.61e-02*(t(l, i)-tm)+8.9e-01))
!                    do not allow ice above 0 celsius unless the top
!                    of the cloud is glaciated and force homogeneous
!                    nucleation at -40 celsius.
                     if ( (t(l, i) >  tm).and.                          &
     &                    (.not.l_glaciated_top(l)) ) then
                       liq_frac_conv(l)=1.0e+00
                     else if (t(l, i) <  tm-4.0e+01) then
                       liq_frac_conv(l)=0.0e+00
                     endif
                  endif
                  condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)      &
     &               =cccwp(lg)*liq_frac_conv(l)                        &
     &               /(total_mass(l)+tiny(cccwp))
                  condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)        &
     &               =cccwp(lg)*(1.0e+00-liq_frac_conv(l))              &
     &               /(total_mass(l)+tiny(cccwp))
                  condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)      &
     &               =cccwp(lg)*liq_frac_conv(l)                        &
     &               /(total_mass(l)+tiny(cccwp))
                  condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)        &
     &               =cccwp(lg)*(1.0e+00-liq_frac_conv(l))              &
     &               /(total_mass(l)+tiny(cccwp))
               endif
            enddo
         enddo
!
!
!        stratiform clouds:
!
!        partition between ice and water depending on the
!        local temperature.
!
         do i=1, nclds
            if (.not.l_cloud_water_partition)                           &
! DEPENDS ON: lsp_focwwil
     &        call lsp_focwwil(t(1, n_layer+1-i), n_profile, liq_frac)
            do l=1, n_profile
               lg=i_gather(l)
               if (l_cloud_water_partition) then
                 if (lca_area(lg, i) > epsilon(lca_area)) then
                   liq_frac(l) = lccwc1(lg, i)                          &
     &               / (lccwc1(lg, i)+lccwc2(lg, i))
                 else
                   liq_frac(l) = 0.0
                 endif
               endif
               if (lca_area(lg, i) >  epsilon(lca_area)) then
                 condensed_mix_rat_area(l, n_layer+1-i                  &
     &              , ip_clcmp_st_water)                                &
     &               =(lccwc1(lg, i)+lccwc2(lg, i))                     &
     &               *liq_frac(l)/lca_area(lg, i)
                 condensed_mix_rat_area(l, n_layer+1-i                  &
     &              , ip_clcmp_st_ice)                                  &
     &               =(lccwc1(lg, i)+lccwc2(lg, i))                     &
     &               *(1.0e+00-liq_frac(l))/lca_area(lg, i)
               else
                 condensed_mix_rat_area(l, n_layer+1-i                  &
     &              , ip_clcmp_st_water)=0.0e+00
                 condensed_mix_rat_area(l, n_layer+1-i                  &
     &              , ip_clcmp_st_ice)=0.0e+00
               endif
!
               if (lca_bulk(lg, i) >  epsilon(lca_bulk)) then
                 condensed_mix_rat_bulk(l, n_layer+1-i                  &
     &              , ip_clcmp_st_water)                                &
     &               =(lccwc1(lg, i)+lccwc2(lg, i))                     &
     &               *liq_frac(l)/lca_bulk(lg, i)
                 condensed_mix_rat_bulk(l, n_layer+1-i                  &
     &              , ip_clcmp_st_ice)                                  &
     &               =(lccwc1(lg, i)+lccwc2(lg, i))                     &
     &               *(1.0e+00-liq_frac(l))/lca_bulk(lg, i)
               else
                 condensed_mix_rat_bulk(l, n_layer+1-i                  &
     &              , ip_clcmp_st_water)=0.0e+00
                 condensed_mix_rat_bulk(l, n_layer+1-i                  &
     &              , ip_clcmp_st_ice)=0.0e+00
               endif
            enddo
         enddo
!
!
!        cloud fractions:
!
       if (lcv_3d_cca) then
         do i=1, nclds
            do l=1, n_profile
               lg=i_gather(l)
               w_cloud(l, n_layer+1-i)                                  &
     &            =cca(lg,i)+(1.0e+00-cca(lg,i))*lca_area(lg, i)
               frac_cloud(l, n_layer+1-i, ip_cloud_type_conv)           &
     &            =cca(lg,i)/(w_cloud(l, n_layer+1-i)+tiny(cca))
               frac_cloud(l, n_layer+1-i, ip_cloud_type_strat)          &
     &            =1.0e+00-frac_cloud(l, n_layer+1-i                    &
     &            , ip_cloud_type_conv)
            enddo
         enddo
       else
         do i=1, nclds
            do l=1, n_profile
              lg=i_gather(l)
               if ( (i <= cct(lg)-1).and.(i >= ccb(lg)) ) then
                  w_cloud(l, n_layer+1-i)                               &
     &               =cca(lg,1)+(1.0e+00-cca(lg,1))*lca_area(lg, i)
                  frac_cloud(l, n_layer+1-i, ip_cloud_type_conv)        &
     &               =cca(lg,1)/(w_cloud(l, n_layer+1-i)+tiny(cca))
               else
                  w_cloud(l, n_layer+1-i)=lca_area(lg, i)
                  frac_cloud(l, n_layer+1-i, ip_cloud_type_conv)        &
     &               =0.0e+00
               endif
               frac_cloud(l, n_layer+1-i, ip_cloud_type_strat)          &
     &            =1.0e+00-frac_cloud(l, n_layer+1-i                    &
     &            , ip_cloud_type_conv)
            enddo
         enddo
       endif
!
!
!
!
      else if (i_cloud_representation == ip_cloud_csiw) then
!
         n_cloud_type=4
!
!        here the clouds are split into four separate types.
!        the partitioning between ice and water is regarded as
!        determining the areas within the grid_box covered by
!        ice or water cloud, rather than as determining the in-cloud
!        mixing ratios. the grid-box mean ice water contents in
!        stratiform clouds may be predicted by the ice microphysics
!        scheme or may be determined as a function of the temperature
!        (lsp_focwwil). in convective clouds the partitioning may be
!        done using the same function, lsp_focwwil, based on a single
!        temperature, or using a partition based on the local
!        temperature.
!
!        convective cloud:
!
          do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)=0.0e+00
               condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)=0.0e+00
               condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)=0.0e+00
               condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)=0.0e+00
            enddo
         enddo
!
         do l=1, n_profile
            total_mass(l)=0.0e+00
         enddo
!
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               lg=i_gather(l)
               if ( (cct(lg) >= n_layer+2-i).and.                       &
     &              (ccb(lg) <= n_layer+1-i) ) then
                  total_mass(l)=total_mass(l)+d_mass(l, i)
               endif
            enddo
         enddo
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               lg=i_gather(l)
               if ( (cct(lg) >= n_layer+2-i).and.                       &
     &              (ccb(lg) <= n_layer+1-i) ) then
                  condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)      &
     &               =cccwp(lg)/(total_mass(l)+tiny(cccwp))
                  condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)        &
     &               =condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)
                  condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)      &
     &               =cccwp(lg)/(total_mass(l)+tiny(cccwp))
                  condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)        &
     &               =condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)
               endif
            enddo
         enddo
!
!        stratiform clouds:
!
         do i=1, nclds
            do l=1, n_profile
               lg=i_gather(l)
               if (lca_area(lg, i) >  epsilon(lca_area)) then
                 condensed_mix_rat_area(l, n_layer+1-i                  &
     &              , ip_clcmp_st_water)                                &
     &              =(lccwc1(lg, i)+lccwc2(lg, i))/lca_area(lg, i)
                 condensed_mix_rat_area(l, n_layer+1-i                  &
     &              , ip_clcmp_st_ice)                                  &
     &              =condensed_mix_rat_area(l, n_layer+1-i              &
     &              , ip_clcmp_st_water)
               else
                 condensed_mix_rat_area(l, n_layer+1-i                  &
     &              , ip_clcmp_st_water)=0.0e+00
                 condensed_mix_rat_area(l, n_layer+1-i                  &
     &              , ip_clcmp_st_ice)=0.0e+00
               endif
!
               if (lca_bulk(lg, i) >  epsilon(lca_bulk)) then
                 condensed_mix_rat_bulk(l, n_layer+1-i                  &
     &              , ip_clcmp_st_water)                                &
     &              =(lccwc1(lg, i)+lccwc2(lg, i))/lca_bulk(lg, i)
                 condensed_mix_rat_bulk(l, n_layer+1-i                  &
     &              , ip_clcmp_st_ice)                                  &
     &              =condensed_mix_rat_bulk(l, n_layer+1-i              &
     &              , ip_clcmp_st_water)
               else
                 condensed_mix_rat_bulk(l, n_layer+1-i                  &
     &              , ip_clcmp_st_water)=0.0e+00
                 condensed_mix_rat_bulk(l, n_layer+1-i                  &
     &              , ip_clcmp_st_ice)=0.0e+00
               endif
            enddo
         enddo
!
!
!        cloud fractions:
!
         if (l_local_cnv_partition) then
!
!           partition between ice and water using the relationships
!           given in bower et al. (1996, q.j. 122 p 1815-1844). ice
!           is allowed in a layer warmer than the freezing point
!           only if the top of the cloud is glaciated.
!
            do l=1, n_profile
!              min is required since cct may be 0 if there is no
!              convective cloud.
               l_glaciated_top(l)                                       &
     &            =(t(l, min(n_layer+2-cct(i_gather(l))                 &
     &            , n_layer-nclds+1)) <  tm)
            enddo

         else
!
!           partition between ice and water as directed by the
!           temperature in the middle of the top layer of the cloud.
!           the partitioning may be precalculated in this case.
!
            do l=1, n_profile
               t_gather(l)=t(l, min(n_layer+2-cct(i_gather(l))          &
     &            , n_layer-nclds+1))
            enddo
! DEPENDS ON: lsp_focwwil
            call lsp_focwwil(t_gather, n_profile, liq_frac_conv)
!
         endif
!
!
         do i=n_layer+1-nclds, n_layer
!
            if (.not. l_cloud_water_partition)                          &
!           partition stratiform clouds using the local temperature.
! DEPENDS ON: lsp_focwwil
     &        call lsp_focwwil(t(1, i), n_profile, liq_frac)
!
          if (lcv_3d_cca) then
            do l=1, n_profile
               lg=i_gather(l)
              convective_cloud_layer(l)=cca(lg,n_layer+1-i)
            enddo
          else
            do l=1, n_profile
            lg=i_gather(l)
               if ( (cct(lg) >= n_layer+2-i).and.                       &
     &              (ccb(lg) <= n_layer+1-i) ) then
                convective_cloud_layer(l)=cca(lg,1)
               else
                  convective_cloud_layer(l)=0.0e+00
               endif
            enddo
          endif
!
            do l=1, n_profile
            lg=i_gather(l)
               w_cloud(l, i)                                            &
     &            =convective_cloud_layer(l)                            &
     &            +(1.0e+00-convective_cloud_layer(l))                  &
     &            *lca_area(lg, n_layer+1-i)
!
               if (l_cloud_water_partition) then
!  partition stratiform clouds using the ratio of cloud water contents.
                 if (lca_area(lg, n_layer+1-i) >                        &
     &             epsilon(lca_area)) then
                   liq_frac(l) = lccwc1(lg, n_layer+1-i) /              &
     &              (lccwc1(lg, n_layer+1-i) + lccwc2(lg, n_layer+1-i))
                 else
                   liq_frac(l) = 0.0e+00
                 endif
               endif
!
               if (l_local_cnv_partition) then
!
!                the partitioning between ice and water must be
!                recalculated for this layer as a function of the
!                local temperature, but ice is allowed above the
!                freezing point only if the top of the cloud is
!                glaciated.
                 liq_frac_conv(l)=max(0.0e+00, min(1.0e+00              &
     &              , 1.61e-02*(t(l, i)-tm)+8.9e-01))
!                do not allow ice above 0 celsius unless the top
!                of the cloud is glaciated and force homogeneous
!                nucleation at -40 celsius.
                 if ( (t(l, i) >  tm).and.                              &
     &                (.not.l_glaciated_top(l)) ) then
                    liq_frac_conv(l)=1.0e+00
                 else if (t(l, i) <  tm-4.0e+01) then
                    liq_frac_conv=0.0e+00
                 endif

               endif
!
               frac_cloud(l, i, ip_cloud_type_sw)                       &
     &            =liq_frac(l)*(1.0e+00-convective_cloud_layer(l))      &
     &            *lca_area(lg, n_layer+1-i)                            &
     &            /(w_cloud(l, i)+tiny(w_cloud))
               frac_cloud(l, i, ip_cloud_type_si)                       &
     &            =(1.0e+00-liq_frac(l))                                &
     &            *(1.0e+00-convective_cloud_layer(l))                  &
     &            *lca_area(lg, n_layer+1-i)                            &
     &            /(w_cloud(l, i)+tiny(w_cloud))
               frac_cloud(l, i, ip_cloud_type_cw)                       &
     &            =liq_frac_conv(l)*convective_cloud_layer(l)           &
     &            /(w_cloud(l, i)+tiny(w_cloud))
               frac_cloud(l, i, ip_cloud_type_ci)                       &
     &            =(1.0e+00-liq_frac_conv(l))*convective_cloud_layer(l) &
     &            /(w_cloud(l, i)+tiny(w_cloud))
!
            enddo
         enddo
!
!
      endif
!
!
!
!     effective radii of water clouds: a microphysical parametrization
!     is available; otherwise standard values are used.
!
      if (l_microphysics) then
!
!        standard values are used for ice crystals, but
!        a parametrization provided by umist and mrf
!        is used for droplets.
!
!        calculate the density of air.
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               density_air(l, i)=p(l, i)/(r*t(l, i))
            enddo
         enddo
!
         do l=1, n_profile
            cc_depth(l)=0.0e+00
         enddo
!
         do l=1, n_profile
            lg=i_gather(l)
!           this loop should be safe even when convective
!           cloud is not present, since ccb should not exceed cct.
            do i=n_layer+2-cct(lg), n_layer+1-ccb(lg)
               cc_depth(l)=cc_depth(l)+d_mass(l, i)/density_air(l, i)
            enddo
         enddo
!
         do l=1, n_profile
            lying_snow_g(l)=lying_snow(i_gather(l))
         enddo
!
! DEPENDS ON: r2_re_mrf_umist
         call r2_re_mrf_umist(n_profile, n_layer, nclds                 &
     &      , i_gather                                                  &
     &      , l_aerosol_ccn                                             &
     &      , sea_salt_film, sea_salt_jet                               &
     &      , l_seasalt_ccn, salt_dim_a, salt_dim_b                     &
     &      , accum_sulphate, diss_sulphate                             &
     &      , aitken_sulphate, lying_snow_g                             &
     &      , i_cloud_representation                                    &
     &      , land_g, flandg_g,  density_air                            &
     &      , condensed_mix_rat_bulk, cc_depth                          &
     &      , condensed_dim_char                                        &
     &      , ntot_diag_g, Ntot_land, Ntot_sea                          &
     &      , strat_lwc_diag_g                                          &
     &      , so4_ccn_diag_g                                            &
     &      , sulp_dim1, sulp_dim2                                      &
     &      , nd_field, nd_profile, nd_layer, id_ct, nd_aerosol_species &
     &      , nd_cloud_component                                        &
     &      )
!
!        constrain the sizes of droplets to lie within the range of
!        validity of the parametrization scheme.
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               condensed_dim_char(l, i, ip_clcmp_st_water)              &
     &            =max(condensed_min_dim(ip_clcmp_st_water)             &
     &            , min(condensed_max_dim(ip_clcmp_st_water)            &
     &            , condensed_dim_char(l, i, ip_clcmp_st_water)))
               condensed_dim_char(l, i, ip_clcmp_cnv_water)             &
     &            =max(condensed_min_dim(ip_clcmp_cnv_water)            &
     &            , min(condensed_max_dim(ip_clcmp_cnv_water)           &
     &            , condensed_dim_char(l, i, ip_clcmp_cnv_water)))
            enddo
         enddo
!
!
!        set microphysical diagnostics. weights for cloud calculated
!        here are used solely for the microphysics and do not have
!        an independent meaning.
!
         if (wgt_conv_flag) then
            if (i_cloud_representation == ip_cloud_conv_strat) then
               do i=1, nclds
                  do l=1, n_profile
                     wgt_conv(col_list(l), row_list(l), i)              &
     &                  =w_cloud(l, n_layer+1-i)                        &
     &                  *frac_cloud(l, n_layer+1-i, ip_cloud_type_conv)
                  enddo
               enddo
            else if (i_cloud_representation == ip_cloud_csiw) then
               do i=1, nclds
                  do l=1, n_profile
                     wgt_conv(col_list(l), row_list(l), i)              &
     &                  =w_cloud(l, n_layer+1-i)                        &
     &                  *frac_cloud(l, n_layer+1-i, ip_cloud_type_cw)
                  enddo
               enddo
            endif
         endif
!
         if (re_conv_flag) then
            do i=1, nclds
               do l=1, n_profile
!                 effective radii are given in microns.
                  re_conv(col_list(l), row_list(l), i)                  &
     &               =condensed_dim_char(l, n_layer+1-i                 &
     &               , ip_clcmp_cnv_water)                              &
     &               *wgt_conv(col_list(l), row_list(l), i)*1.0e+06
               enddo
            enddo
         endif
!
         if (wgt_strat_flag) then
            if (i_cloud_representation == ip_cloud_conv_strat) then
               do i=1, nclds
                  do l=1, n_profile
                     wgt_strat(col_list(l), row_list(l), i)             &
     &                  =w_cloud(l, n_layer+1-i)                        &
     &                  *frac_cloud(l, n_layer+1-i                      &
     &                  , ip_cloud_type_strat)
                  enddo
               enddo
            else if (i_cloud_representation == ip_cloud_csiw) then
               do i=1, nclds
                  do l=1, n_profile
                     wgt_strat(col_list(l), row_list(l), i)             &
     &                  =w_cloud(l, n_layer+1-i)                        &
     &                  *frac_cloud(l, n_layer+1-i, ip_cloud_type_sw)
                  enddo
               enddo
            endif
         endif
!
         if (re_strat_flag) then
            do i=1, nclds
               do l=1, n_profile
!                 effective radii are given in microns.
                  re_strat(col_list(l), row_list(l), i)                 &
     &               =condensed_dim_char(l, n_layer+1-i                 &
     &               , ip_clcmp_st_water)                               &
     &               *wgt_strat(col_list(l), row_list(l), i)*1.0e+06
               enddo
            enddo
         endif

         if (lwp_strat_flag) then
            do i=1, nclds
               do l=1, n_profile
                  lwp_strat(col_list(l), row_list(l), i)                &
     &               =condensed_mix_rat_area(l, n_layer+1-i             &
     &               , ip_clcmp_st_water)*d_mass(l, n_layer+1-i)        &
     &               *wgt_strat(col_list(l), row_list(l), i)
               enddo
            enddo
         endif

         if (ntot_diag_flag) then
            do i=1, nclds
               do l=1, n_profile
                  ntot_diag(col_list(l), row_list(l), i)                &
     &               =ntot_diag_g(l, n_layer+1-i)                       &
     &               *wgt_strat(col_list(l), row_list(l), i)
               enddo
            enddo
         endif

         if (strat_lwc_diag_flag) then
            do i=1, nclds
               do l=1, n_profile
                  strat_lwc_diag(col_list(l), row_list(l), i)           &
     &               =strat_lwc_diag_g(l, n_layer+1-i)                  &
     &               *wgt_strat(col_list(l), row_list(l), i)
               enddo
            enddo
         endif

! non-cloud diagnostics are "weighted" by the conditional sampling
! weight cond_samp_wgt, but as this is 1.0 if the sw radiation is
! active, and 0.0 if it is not, there is no need to actually
! multiply by it.

         if (cond_samp_wgt_flag) then
            do i=1, nclds
               do l=1, n_profile
                  cond_samp_wgt(col_list(l), row_list(l), i)=1.0
               enddo
            enddo
         endif

         if (so4_ccn_diag_flag) then
            do i=1, nclds
               do l=1, n_profile
                  so4_ccn_diag(col_list(l), row_list(l), i)             &
     &                    =so4_ccn_diag_g(l, n_layer+1-i)
               enddo
            enddo
         endif
!
!
      else
!
!        all effective radii are set to standard values.
!
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               condensed_dim_char(l, i, ip_clcmp_st_water)=7.e-6
               condensed_dim_char(l, i, ip_clcmp_cnv_water)=7.e-6
            enddo
         enddo
!
      endif
!
!
!
!     set the characteristic dimensions of ice crystals:
!
!     ice crystals in stratiform clouds:
!
      if (i_condensed_param(ip_clcmp_st_ice) ==                         &
     &   ip_slingo_schrecker_ice) then
!
!        this parametrization is based on the effective radius
!        and a standard value of 30-microns is assumed.
!
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               condensed_dim_char(l, i, ip_clcmp_st_ice)=30.e-6
            enddo
         enddo
!
      else if (                                                         &
     &  (i_condensed_param(ip_clcmp_st_ice) == ip_ice_adt) .OR.         &
     &  (i_condensed_param(ip_clcmp_st_ice) == ip_ice_adt_10)           &
     &  ) then
!
!        this parametrization is based on the mean maximum
!        dimension of the crystal, determined as a function of
!        the local temperature. the size is limited to its value
!        at the freezing level.
!
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               condensed_dim_char(l, i, ip_clcmp_st_ice)                &
     &            =min(7.198755e-04                                     &
     &            , exp(5.522e-02*(t(l, i)-2.7965e+02))/9.702e+02)
            enddo
         enddo
      else if (i_condensed_param(IP_clcmp_st_ice) ==                    &
     &   IP_ice_fu_phf) then
!
!      Aggregate parametrization based on effective dimension.
!      In the initial form, the same approach is used for stratiform
!      and convective cloud.
!
!      The fit provided here is based on Stephan Havemann's fit of
!      Dge with temperature, consistent with David Mitchell's treatment
!      of the variation of the size distribution with temperature. The
!      parametrization of the optical properties is based on De
!      (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
!      (=(2*SQRT(3)/3)*volume/projected area), which explains the
!      conversion factor. The fit to Dge is in two sections, because
!      Mitchell's relationship predicts a cusp at 216.208 K. Limits
!      of 8 and 124 microns are imposed on Dge: these are based on this
!      relationship and should be reviewed if it is changed. Note also
!      that the relationship given here is for polycrystals only.
       do i=n_layer+1-nclds, n_layer
         do l=1, n_profile
!          Preliminary calculation of Dge.
           IF (t(l, i) < t_switch) THEN
             condensed_dim_char(l, i, IP_clcmp_st_ice)                  &
     &         = a0_agg_cold*EXP(s0_agg*(t(l, i)-t0_agg))+b0_agg_cold
           ELSE
             condensed_dim_char(l, i, IP_clcmp_st_ice)                  &
     &         = a0_agg_warm*EXP(s0_agg*(t(l, i)-t0_agg))+b0_agg_warm
           ENDIF
!          Limit and convert to De.
           condensed_dim_char(l, i, IP_clcmp_st_ice)                    &
     &       = (3.0/2.0)*(3.0/(2.0*SQRT(3.0)))*                         &
     &         MIN(1.24E-04, MAX(8.0E-06,                               &
     &         condensed_dim_char(l, i, IP_clcmp_st_ice)))
         enddo
       enddo
!
      endif
!
!
!     ice crystals in convective clouds:
!
      if (i_condensed_param(ip_clcmp_cnv_ice) ==                        &
     &   ip_slingo_schrecker_ice) then
!
!        this parametrization is based on the effective radius
!        and a standard value of 30-microns is assumed.
!
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               condensed_dim_char(l, i, ip_clcmp_cnv_ice)=30.e-6
            enddo
         enddo
!
      else if (                                                         &
     &  (i_condensed_param(ip_clcmp_cnv_ice) == ip_ice_adt) .OR.        &
     &  (i_condensed_param(ip_clcmp_cnv_ice) == ip_ice_adt_10)          &
     &  ) then
!
!        this parametrization is based on the mean maximum
!        dimension of the crystal, determined as a function of
!        the local temperature. the size is limited to its value
!        at the freezing level.
!
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               condensed_dim_char(l, i, ip_clcmp_cnv_ice)               &
     &            =min(7.198755e-04                                     &
     &            , exp(5.522e-02*(t(l, i)-2.7965e+02))/9.702e+02)
            enddo
         enddo
!
      else if (i_condensed_param(IP_clcmp_cnv_ice) ==                   &
     &   IP_ice_fu_phf) then
!
!      Aggregate parametrization based on effective dimension.
!      In the initial form, the same approach is used for stratiform
!      and convective cloud.
!
!      The fit provided here is based on Stephan Havemann's fit of
!      Dge with temperature, consistent with David Mitchell's treatment
!      of the variation of the size distribution with temperature. The
!      parametrization of the optical properties is based on De
!      (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
!      (=(2*SQRT(3)/3)*volume/projected area), which explains the
!      conversion factor. The fit to Dge is in two sections, because
!      Mitchell's relationship predicts a cusp at 216.208 K. Limits
!      of 8 and 124 microns are imposed on Dge: these are based on this
!      relationship and should be reviewed if it is changed. Note also
!      that the relationship given here is for polycrystals only.
       do i=n_layer+1-nclds, n_layer
         do l=1, n_profile
!          Preliminary calculation of Dge.
           if (t(l, i) < t_switch) then
             condensed_dim_char(l, i, IP_clcmp_cnv_ice)                 &
     &         = a0_agg_cold*EXP(s0_agg*(t(l, i)-t0_agg))+b0_agg_cold
           else
             condensed_dim_char(l, i, IP_clcmp_cnv_ice)                 &
     &         = a0_agg_warm*EXP(s0_agg*(t(l, i)-t0_agg))+b0_agg_warm
           endif
!          Limit and convert to De.
           condensed_dim_char(l, i, IP_clcmp_cnv_ice)                   &
     &       = (3.0/2.0)*(3.0/(2.0*SQRT(3.0)))*                         &
     &         MIN(1.24E-04, MAX(8.0E-06,                               &
     &         condensed_dim_char(l, i, IP_clcmp_cnv_ice)))
         enddo
       enddo
!
      endif
!
!
!
!     constrain the sizes of ice crystals to lie within the range
!     of validity of the parametrization scheme.
      do i=n_layer+1-nclds, n_layer
         do l=1, n_profile
            condensed_dim_char(l, i, ip_clcmp_st_ice)                   &
     &         =max(condensed_min_dim(ip_clcmp_st_ice)                  &
     &         , min(condensed_max_dim(ip_clcmp_st_ice)                 &
     &         , condensed_dim_char(l, i, ip_clcmp_st_ice)))
            condensed_dim_char(l, i, ip_clcmp_cnv_ice)                  &
     &         =max(condensed_min_dim(ip_clcmp_cnv_ice)                 &
     &         , min(condensed_max_dim(ip_clcmp_cnv_ice)                &
     &         , condensed_dim_char(l, i, ip_clcmp_cnv_ice)))
         enddo
      enddo
!
!
!
      return
      END SUBROUTINE r2_set_cloud_field
#endif
#endif
