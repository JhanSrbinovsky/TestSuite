#if defined(A01_3Z)
! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************
!+ Subroutine to calculate an observed effective radius.
!
! Purpose:
!   an effective radius as observed from above the cloud-top is
!   calculated.
!
! Method:
!  For each type of cloud containing water in any layer the effective
!  radius is weighted with the product of the area of the cloud and the
!  probability that light emitted from the cloud reaches the observing
!  instrument.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             13-02-06                New version
!                                               (J.-C. Thelen)
!
! Description of code:
!   Fortran 90 in fixed format
!
!---------------------------------------------------------------------
      subroutine r2_cloud_level_diag(ierr, n_profile, n_layer, nclds    &
     &   , i_gather                                                     &
     &   , i_cloud, i_cloud_representation                              &
     &   , w_cloud, frac_cloud                                          &
     &   , condensed_mix_ratio, condensed_re                            &
     &   , l_observed_re, weighted_re, sum_weight_re                    &
     &   , col_list, row_list, row_length, rows                         &
     &   , nd_field, nd_profile, id_ct, nd_layer                        &
     &   , nd_cloud_component, nd_cloud_type                            &
     &   )
!
!
!
      implicit none
!
!
!     comdecks included.
#include "error_pcf3z.h"
#include "def_std_io_icf3z.h"
#include "cloud_component_pcf3z.h"
#include "cloud_type_pcf3z.h"
#include "cloud_representation_pcf3z.h"
#include "cloud_scheme_pcf3z.h"
!
!
      integer                                                           &
                !, intent(out)
     &     ierr
!             error flag
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
!             size of array of arrays passed from main code
     &   , nd_profile                                                   &
!             size of array of profiles
     &   , id_ct                                                        &
!             Topmost declared cloudy layer
     &   , nd_layer                                                     &
!             maximum number of layers
     &   , nd_cloud_component                                           &
!            Size allocated for condensed components
     &   , nd_cloud_type
!            Size allocated for number of types of cloud
!
!     actual sizes used:
      integer                                                           &
                !, intent(in)
     &     n_profile                                                    &
!             number of profiles
     &   , n_layer                                                      &
!             number of atmospheric layers in radiation
     &   , nclds                                                        &
!             number of cloudy levels
     &   , i_gather(nd_field)
!             list of gathered points
      integer, intent(in) :: col_list(nd_field)
!                              ew indices of gathered points in the
!                              2-d domain
      integer, intent(in) :: row_list(nd_field)
!                              ns indices of gathered points in the
!                              2-d domain
!
!     logical flags for diagnostics
      logical                                                           &
                !, intent(in)
     &     l_observed_re
!             flag to enable diagnosis of effective radius seen from
!             space (n.b. the routine is at present called only if
!             this is true, but its presence here allows for possible
!             future extension of the routine).
!
!     representation of clouds
      integer                                                           &
                !, intent(in)
     &     i_cloud_representation                                       &
!             representation of clouds
     &   , i_cloud
!             treatment of overlaps
!
      real                                                              &
                !, intent(in)
     &     w_cloud(nd_profile, nd_layer)                                &
!             total amounts of cloud
     &   , frac_cloud(nd_profile, nd_layer, nd_cloud_type)              &
!             fraction of types of cloud
     &   , condensed_re(nd_profile,id_ct:nd_layer,nd_cloud_component)   &
!             effective radii of cloudy components
     &   , condensed_mix_ratio(nd_profile, id_ct: nd_layer              &
     &      , nd_cloud_component)
!             mass mixing ratios of condensed components
!
      real                                                              &
                !, intent(out)
     &     weighted_re(row_length, rows)                                &
!             weighted sum of effective radius and weighting function
     &   , sum_weight_re(row_length, rows)
!             sum of weights for effective radius
!
!
!
!     local variables:
      integer                                                           &
     &     i                                                            &
!             loop variable
     &   , l                                                            &
!             loop variable
     &   , i_inv
!             inverted loop index
      real                                                              &
     &     trans_overlying_space(nd_profile)                            &
!             probability of a photon in clear air in the level above
!             the current one reaching space
     &   , area_exposed(nd_profile)                                     &
!             total area of cloud in the current layer exposed to
!             clear air in the layer above
     &   , area_exposed_st(nd_profile)                                  &
!             total area of stratiform cloud in the current layer
!             exposed to clear air in the layer above
     &   , area_exposed_cnv(nd_profile)                                 &
!             total area of convective cloud in the current layer
!             exposed to clear air in the layer above
     &   , area_clear_above(nd_profile)                                 &
!             area of the clear sky region in the layer above
     &   , area_strat(nd_profile)                                       &
!             area of stratiform cloud in the current layer
     &   , area_strat_above(nd_profile)                                 &
!             area of stratiform cloud in the layer above
     &   , area_conv(nd_profile)                                        &
!             area of convective cloud in the current layer
     &   , area_conv_above(nd_profile)                                  &
!             area of convective cloud in the layer above
     &   , area_clear_clear(nd_profile)                                 &
!             area of boundary where clear sky overlies clear sky
     &   , area_clear(nd_profile)                                       &
!             area of clear sky in the current layer
!             down to a level
     &   , area_uncorrelated(nd_profile)                                &
!             uncorrelated region on the interface
     &   , weighted_re_g(nd_profile)                                    &
!             weighted sum of effective radius and weighting function
     &   , sum_weight_re_g(nd_profile)
!             sum of weights for effective radius
!
!     variables for gathering
      integer                                                           &
     &     n_list                                                       &
!             number of points in list
     &   , l_list(nd_profile)
!             indices of points in list
!
!     indicator function
      real                                                              &
     &     chi_cnv(nd_profile)                                          &
!             convective indicator function
     &   , chi_st(nd_profile)
!             stratiform indicator function
!
!
!
!
!     initialization of local fields.
!
      if (l_observed_re) then
         do l=1, n_profile
            weighted_re_g(l)=0.0e+00
            sum_weight_re_g(l)=0.0e+00
         enddo
      endif
!
!     initialize the transmision above clouds.
      do l=1, n_profile
         trans_overlying_space(l)=1.0e+00
         area_clear_above(l)=1.0e+00
      enddo
      if (l_observed_re.and.(i_cloud == ip_cloud_triple)) then
         do l=1, n_profile
            area_strat_above(l)=0.0e+00
            area_conv_above(l)=0.0e+00
         enddo
      endif
!
!     step down through the atmosphere calculating contributions to
!     the diagnostics and subsequently allowing for transmission
!     through the current layer.
!
      do i=nclds, 1, -1
         i_inv=n_layer+1-i
!
         do l=1, n_profile
            area_clear(l)=1.0e+00-w_cloud(l, i_inv)
         enddo
!
!        calculate the local area of cloud radiating into clear air.
         if (i_cloud == ip_cloud_mix_random) then
            do l=1, n_profile
               area_exposed(l)=w_cloud(l, i_inv)                        &
     &            *area_clear_above(l)
            enddo
         else if ( (i_cloud == ip_cloud_mix_max).or.                    &
     &             (i_cloud == ip_cloud_triple) ) then
            do l=1, n_profile
               area_exposed(l)=max(0.0e+00, (w_cloud(l, i_inv)          &
     &            +area_clear_above(l)-1.0e+00))
            enddo
         endif
!
!
!
!
         if (l_observed_re) then
!
!
!
            if (i_cloud_representation == ip_cloud_conv_strat) then
!
!
               if ( (i_cloud == ip_cloud_mix_max).or.                   &
     &              (i_cloud == ip_cloud_mix_random) ) then
!
!                 if the overlap of convective cloud is not assumed
!                 to be coherent the overall exposed area may be
!                 partitioned according to the fractional
!                 contributions of cloud in the current layer.
!
                  do l=1, n_profile
                     area_exposed_st(l)=area_exposed(l)                 &
     &                  *frac_cloud(l, i_inv, ip_cloud_type_strat)
                     area_exposed_cnv(l)=area_exposed(l)                &
     &                  *frac_cloud(l, i_inv, ip_cloud_type_conv)
                  enddo
!
               else if (i_cloud == ip_cloud_triple) then
!
!                 here, the different types of clouds overlap
!                 coherently so stratiform cloud will be exposed
!                 only if there is less stratiform cloud in the
!                 layer above and more clear air in the layer above:
!                 under these conditions the non-correlated areas
!                 overlap randomly.
!
                  do l=1, n_profile
                     area_strat(l)=w_cloud(l, i_inv)                    &
     &                  *frac_cloud(l, i_inv, ip_cloud_type_strat)
                     area_conv(l)=w_cloud(l, i_inv)                     &
     &                 *frac_cloud(l, i_inv, ip_cloud_type_conv)
                     area_uncorrelated(l)=1.0e+00                       &
     &                 -min(area_clear(l), area_clear_above(l))         &
     &                 -min(area_strat(l), area_strat_above(l))         &
     &                 -min(area_conv(l), area_conv_above(l))
!                    first find the area of uncorrelated
!                    stratiform cloud.
                     area_exposed_st(l)=max(0.0e+00                     &
     &                  , (area_strat(l)-area_strat_above(l)))
                     area_exposed_st(l)=max(0.0e+00, area_exposed_st(l) &
     &                  *(area_clear_above(l)-area_clear(l)))
!                    now normalize within the uncorrelated region.
!                    if the uncorrelated area is 0 the exposed area
!                    must be 0, so no second branch of the if-test
!                    is required.
                     if (area_uncorrelated(l) >  0.0e+00)               &
     &                 area_exposed_st(l)                               &
     &                   =area_exposed_st(l)/area_uncorrelated(l)
                     area_exposed_cnv(l)                                &
     &                  =area_exposed(l)-area_exposed_st(l)
                  enddo
               else
                  write(iu_err, '(/a)')                                 &
     &              '*** error: the diagnostic of observed re has not ' &
     &              //'been implemented with this overlap option.'
                  ierr=i_err_fatal
                  return
               endif
!
!              the indicator functions for liquid water in
!              convective or straiform clouds are set to 1
!              if there is any liquid water and to 0 otherwise.
               do l=1, n_profile
                  if (condensed_mix_ratio(l, i_inv, ip_clcmp_cnv_water) &
     &                >  0.0e+00) then
                     chi_cnv(l)=1.0e+00
                  else
                     chi_cnv(l)=0.0e+00
                  endif
                  if (condensed_mix_ratio(l, i_inv, ip_clcmp_st_water)  &
     &                >  0.0e+00) then
                     chi_st(l)=1.0e+00
                  else
                     chi_st(l)=0.0e+00
                  endif
               enddo
!
!              include contributions from convective and stratiform
!              water clouds.
               do l=1, n_profile
                  weighted_re_g(l)=weighted_re_g(l)                     &
     &               +trans_overlying_space(l)                          &
     &               *(area_exposed_cnv(l)*chi_cnv(l)                   &
     &               *condensed_re(l, i_inv, ip_clcmp_cnv_water)        &
     &               +area_exposed_st(l)*chi_st(l)                      &
     &               *condensed_re(l, i_inv, ip_clcmp_st_water))
                  sum_weight_re_g(l)=sum_weight_re_g(l)                 &
     &               +trans_overlying_space(l)                          &
     &               *(area_exposed_cnv(l)*chi_cnv(l)                   &
     &               +area_exposed_st(l)*chi_st(l))
               enddo
!
            else if (i_cloud_representation == ip_cloud_csiw) then
!
               if ( (i_cloud == ip_cloud_mix_max).or.                   &
     &              (i_cloud == ip_cloud_mix_random) ) then
!
!                 if the overlap of convective cloud is not assumed
!                 to be coherent the overall exposed area may be
!                 partitioned according to the fractional
!                 contributions of cloud in the current layer.
!                 the exposed areas include only the parts of the
!                 clouds containing water droplets.
!
                  do l=1, n_profile
                     area_exposed_st(l)=area_exposed(l)                 &
     &                  *frac_cloud(l, i_inv, ip_cloud_type_sw)
                     area_exposed_cnv(l)=area_exposed(l)                &
     &                  *frac_cloud(l, i_inv, ip_cloud_type_cw)
                  enddo
!
               else if (i_cloud == ip_cloud_triple) then
!
!                 here, the different types of clouds overlap
!                 coherently so stratiform cloud will be exposed
!                 only if there is less stratiform cloud in the
!                 layer above and more clear air in the layer above:
!                 under these conditions the non-correlated areas
!                 overlap randomly.
!                 the actual exposed areas of convective or
!                 stratiform cloud must then be weighted by factors
!                 representing the liquid portion of each cloud, since
!                 nothing is retrieved over ice. (the horizontal
!                 arrangement of ice and water within either type of
!                 cloud is random).
!
                  do l=1, n_profile
!
                     area_strat(l)=w_cloud(l, i_inv)                    &
     &                  *(frac_cloud(l, i_inv, ip_cloud_type_sw)        &
     &                  +frac_cloud(l, i_inv, ip_cloud_type_si))
                     area_conv(l)=w_cloud(l, i_inv)                     &
     &                  *(frac_cloud(l, i_inv, ip_cloud_type_cw)        &
     &                  +frac_cloud(l, i_inv, ip_cloud_type_ci))
                     area_uncorrelated(l)=1.0e+00                       &
     &                  -min(area_clear(l), area_clear_above(l))        &
     &                  -min(area_strat(l), area_strat_above(l))        &
     &                  -min(area_conv(l), area_conv_above(l))
                     area_exposed_st(l)=max(0.0e+00                     &
     &                  , (area_strat(l)-area_strat_above(l)))
                     if (area_uncorrelated(l) >  0.0e+00) then
                        area_exposed_st(l)                              &
     &                     =max(0.0e+00, area_exposed_st(l)             &
     &                     *(area_clear_above(l)-area_clear(l)))        &
     &                     /area_uncorrelated(l)
                     else
                        area_exposed_st(l)=0.0e+00
                     endif
                     area_exposed_cnv(l)                                &
     &                  =area_exposed(l)-area_exposed_st(l)
!
                     if (frac_cloud(l, i_inv, ip_cloud_type_cw)         &
     &                   >  0.0e+00) then
                        area_exposed_cnv(l)=area_exposed_cnv(l)         &
     &                     /(1.0e+00                                    &
     &                     +frac_cloud(l, i_inv, ip_cloud_type_ci)      &
     &                     /frac_cloud(l, i_inv, ip_cloud_type_cw))
                     else
                        area_exposed_cnv(l)=0.0e+00
                     endif
!
                     if (frac_cloud(l, i_inv, ip_cloud_type_sw)         &
     &                   >  0.0e+00) then
                        area_exposed_st(l)=area_exposed_st(l)           &
     &                     /(1.0e+00                                    &
     &                     +frac_cloud(l, i_inv, ip_cloud_type_si)      &
     &                     /frac_cloud(l, i_inv, ip_cloud_type_sw))
                     else
                        area_exposed_st(l)=0.0e+00
                     endif
!
                  enddo
               else
                  write(iu_err, '(/a)')                                 &
     &              '*** error: the diagnostic of observed re has not ' &
     &              //'been implemented with this overlap option.'
                  ierr=i_err_fatal
                  return
               endif
!
!
               do l=1, n_profile
!
                  weighted_re_g(l)=weighted_re_g(l)                     &
     &               +trans_overlying_space(l)                          &
     &               *(area_exposed_cnv(l)                              &
     &               *condensed_re(l, i_inv, ip_clcmp_cnv_water)        &
     &               +area_exposed_st(l)                                &
     &               *condensed_re(l, i_inv, ip_clcmp_st_water))
                  sum_weight_re_g(l)=sum_weight_re_g(l)                 &
     &               +trans_overlying_space(l)                          &
     &               *(area_exposed_cnv(l)+area_exposed_st(l))
               enddo
!
            endif
!
!
         endif
!
!
!
!        advance the stored quantities refferring to overlying layers.
!
!
!      the transmission to space currently holds the probability that
!      a photon travelling upwards in the clear air in the layer above
!      will escape to space without encountering a cloud. to advance
!      this to the current layer it must be multiplied by a factor
!      representing the overlap assumption at the top of the present
!      layer.
!
         if (i_cloud == ip_cloud_mix_random) then
!
            do l=1, n_profile
               trans_overlying_space(l)=trans_overlying_space(l)        &
     &            *area_clear_above(l)
            enddo
!
         else if ( (i_cloud == ip_cloud_mix_max).or.                    &
     &             (i_cloud == ip_cloud_triple) ) then
!
            do l=1, n_profile
               area_clear_clear(l)=min(area_clear(l)                    &
     &            , area_clear_above(l))
               if (area_clear(l) >  0.0e+00) then
                  trans_overlying_space(l)=trans_overlying_space(l)     &
     &               *area_clear_clear(l)/area_clear(l)
               else
                  trans_overlying_space(l)=0.0e+00
               endif
            enddo
!
         endif
!
!        advance the areas of cloud.
         do l=1, n_profile
            area_clear_above(l)=area_clear(l)
         enddo
         if (i_cloud_representation == ip_cloud_conv_strat) then
            do l=1, n_profile
               area_strat_above(l)=w_cloud(l, i_inv)                    &
     &            *frac_cloud(l, i_inv, ip_cloud_type_strat)
            enddo
         else if (i_cloud_representation == ip_cloud_csiw) then
            do l=1, n_profile
               area_strat_above(l)=area_strat(l)
               area_conv_above(l)=area_conv(l)
            enddo
         endif
!
      enddo
!
!
!
      if (l_observed_re) then
!        scatter the diagnostics back to the output arrays and convert
!        to microns (to avoid fields being corrupted by packing).
         do l=1, n_profile
            weighted_re(col_list(l), row_list(l))                       &
     &        =1.0e+06*weighted_re_g(l)
            sum_weight_re(col_list(l), row_list(l))                     &
     &        =sum_weight_re_g(l)
         enddo
      endif
!
!
!
      return
      END SUBROUTINE r2_cloud_level_diag
#endif
