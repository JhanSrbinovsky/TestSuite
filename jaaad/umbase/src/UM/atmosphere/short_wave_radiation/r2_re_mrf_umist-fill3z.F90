#if defined(A70_1Z)
#if defined(A01_3Z) ||  defined(A02_3Z)
! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************
!
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
      subroutine r2_re_mrf_umist(n_profile, n_layer, nclds              &
     &   , i_gather                                                     &
     &   , l_aerosol_ccn                                                &
     &   , sea_salt_film, sea_salt_jet                                  &
     &   , l_seasalt_ccn, salt_dim_a, salt_dim_b                        &
     &   , accum_sulphate, diss_sulphate, aitken_sulphate               &
     &   , lying_snow_g                                                 &
     &   , i_cloud_representation                                       &
     &   , land_g, flandg_g, density_air                                &
     &   , condensed_mix_ratio, cc_depth                                &
     &   , condensed_re                                                 &
     &   , ntot_diag_g, Ntot_land, Ntot_sea                             &
     &   , strat_lwc_diag_g                                             &
     &   , so4_ccn_diag_g                                               &
     &   , sulp_dim1, sulp_dim2                                         &
     &   , nd_field, nd_profile, nd_layer, id_ct, nd_aerosol_species    &
     &   , nd_cloud_component                                           &
     &   )
!
!
!
      implicit none
!
!
!     comdecks included:
#include "c_pi.h"
#include "c_densty.h"
#include "c_micro.h"
#include "cloud_component_pcf3z.h"
#include "cloud_representation_pcf3z.h"
!
!
!     dummy arguments:
!
!     sizes of arrays:
      integer                                                           &
                !, intent(in)
     &     nd_field                                                     &
!             size of input fields to the radiation
     &   , nd_profile                                                   &
!             maximum number of profiles
     &   , nd_layer                                                     &
!             maximum number of layers
     &   , nd_aerosol_species                                           &
!             maximum number of aerosol species
     &   , nd_cloud_component                                           &
!             maximum number of components in clouds
     &   , id_ct                                                        &
!             Topmost declared cloudy layer
     &   , sulp_dim1                                                    &
!             1st dimension of arrays of sulphate
     &   , sulp_dim2                                                    &
!             2nd dimension of arrays of sulphate
     &   , salt_dim_a                                                   &
!             1st dimension of arrays of sea-salt
     &   , salt_dim_b
!             2nd dimension of arrays of sea-salt
!
      integer                                                           &
                !, intent(in)
     &     n_profile                                                    &
!             number of atmospheric profiles
     &   , n_layer                                                      &
!             number of layers seen in radiation
     &   , nclds
!             number of cloudy levels
!
      integer                                                           &
                !, intent(in)
     &     i_gather(nd_field)
!             list of points to be gathered
      logical                                                           &
                !, intent(in)
     &     land_g(nd_profile)
!             gathered mask for land points
      integer                                                           &
                !, intent(in)
     &     i_cloud_representation
!             representation of clouds
!
!     variables for aerosols
      logical                                                           &
                !, intent(in)
     &     l_aerosol_ccn                                                &
!             flag to use aerosols to find ccn.
     &   , l_seasalt_ccn
!             flag to use sea-salt parametrization for ccn
      real                                                              &
                !, intent(in)
     &     accum_sulphate(sulp_dim1, sulp_dim2)                         &
!             mixing ratios of accumulation mode sulphate
     &   , aitken_sulphate(sulp_dim1, sulp_dim2)                        &
!             mixing ratios of aitken-mode sulphate
     &   , diss_sulphate(sulp_dim1, sulp_dim2)                          &
!             mixing ratios of dissolved sulphate
     &   , sea_salt_film(salt_dim_a, salt_dim_b)                        &
!             number concentration of film-mode sea-salt aerosol
     &   , sea_salt_jet(salt_dim_a, salt_dim_b)
!             number concentration of jet-mode sea-salt aerosol
!
      real                                                              &
                !, intent(in)
     &     density_air(nd_profile, nd_layer)
!             density of air
!
      real                                                              &
                !, intent(in)
     &     condensed_mix_ratio(nd_profile, id_ct: nd_layer              &
     &        , nd_cloud_component)                                     &
!             mixing ratios of condensed species
     &   , cc_depth(nd_profile)
!             depth of convective cloud
!
      real                                                              &
                !, intent(in)
     &     Ntot_land                                                    &
!             number of droplets over land / m-3
     &   , Ntot_sea
!             number of droplets over sea / m-3

      real                                                              &
                !, intent(out)
     &     condensed_re(nd_profile, id_ct: nd_layer, nd_cloud_component)
!             effective radii of condensed components of clouds
!
      real                                                              &
                !, intent(out)
     &     ntot_diag_g(nd_profile, nd_layer)                            &
!             diagnostic array for ntot (gathered)
     &   , strat_lwc_diag_g(nd_profile, nd_layer)                       &
!             diagnostic array for stratiform lwc (gathered)
     &   , so4_ccn_diag_g(nd_profile, nd_layer)
!             diagnostic array for so4 ccn mass conc (gathered)
!
!
      real                                                              &
     &     lying_snow_g(nd_profile)                                     &
!             gathered snow depth (>5000m = land ice sheet)
     &   , flandg_g(nd_profile)
!             gathered global land fraction
!
!     local variables:
      integer                                                           &
     &     i                                                            &
!             loop variable
     &   , l                                                            &
!             loop variable
     &   , sulphate_ptr_a                                               &
     &   , sulphate_ptr_b                                               &
!             pointers for sulphate arrays
     &   , seasalt_ptr_a                                                &
     &   , seasalt_ptr_b
!             pointers for sea-salt arrays
!
      real                                                              &
     &     total_mix_ratio_st(nd_profile)                               &
!             total mixing ratio of water substance in stratiform cloud
     &   , total_mix_ratio_cnv(nd_profile)
!             total mixing ratio of water substance in stratiform cloud
!
      real                                                              &
     &     n_drop(nd_profile, nd_layer)                                 &
!             number density of droplets
     &   , kparam
!             ratio of cubes of volume radius to effective radius
!
!     fixed constants of the parametrization:
      real                                                              &
     &     deep_convective_cloud
!             threshold value for deep convective cloud
      parameter(                                                        &
     &     deep_convective_cloud=5.0e+02                                &
     &   )
!
      real                                                              &
     &     dummy_2d(1,1)
!             dummy variable
!
!     functions called:
      real                                                              &
     &     number_droplet
!             function to calculate the number of clouds droplets
      external                                                          &
     &     number_droplet
!
!
!
!
!     calculate the number density of droplets
!
      do i=n_layer+1-nclds, n_layer
         do l=1, n_profile
            if (l_aerosol_ccn) then
               sulphate_ptr_a=i_gather(l)
               sulphate_ptr_b=n_layer+1-i
            else
               sulphate_ptr_a=1
               sulphate_ptr_b=1
            endif
            if (l_seasalt_ccn) then
               seasalt_ptr_a=i_gather(l)
               seasalt_ptr_b=n_layer+1-i
            else
               seasalt_ptr_a=1
               seasalt_ptr_b=1
            endif
! DEPENDS ON: number_droplet
            n_drop(l, i)=number_droplet(l_aerosol_ccn, .true.           &
     &         , aitken_sulphate(sulphate_ptr_a, sulphate_ptr_b)        &
     &         , accum_sulphate(sulphate_ptr_a, sulphate_ptr_b)         &
     &         , diss_sulphate(sulphate_ptr_a, sulphate_ptr_b)          &
     &         , l_seasalt_ccn                                          &
     &         , sea_salt_film(seasalt_ptr_a, seasalt_ptr_b)            &
     &         , sea_salt_jet(seasalt_ptr_a, seasalt_ptr_b)             &
     &         , .false., dummy_2d(1,1), .false.                        &
     &         , dummy_2d(1,1), dummy_2d(1,1)                           &
     &         , .false., dummy_2d(1,1), dummy_2d(1,1)                  &
     &         , density_air(l, i)                                      &
     &         , lying_snow_g(l)                                        &
     &         , flandg_g(l)                                            &
     &         , Ntot_land, Ntot_sea                                    &
     &         )
         enddo
      enddo

!  diagnose so4 aerosol concentrations. mass mixing ratio of ammonium
!  sulphate is converted to microgrammes of the sulphate ion per m3
!  for diagnostic purposes.

      if (l_aerosol_ccn) then
         do i=n_layer+1-nclds, n_layer
            do l=1, n_profile
               so4_ccn_diag_g(l, i)=                                    &
     &                  (aitken_sulphate(i_gather(l), n_layer+1-i)      &
     &                   +accum_sulphate(i_gather(l), n_layer+1-i)      &
     &                    +diss_sulphate(i_gather(l), n_layer+1-i))     &
     &                 * density_air(l, i) * (96./132.) * 1.0e+09
            enddo
         enddo
      endif
!
      do i=n_layer+1-nclds, n_layer
!
!        find the total mixing ratio of water substance in the cloud
!        as implied by the representation.
         if (i_cloud_representation == ip_cloud_conv_strat) then
            do l=1, n_profile
               total_mix_ratio_st(l)                                    &
     &            =condensed_mix_ratio(l, i, ip_clcmp_st_water)         &
     &            +condensed_mix_ratio(l, i, ip_clcmp_st_ice)
               total_mix_ratio_cnv(l)                                   &
     &            =condensed_mix_ratio(l, i, ip_clcmp_cnv_water)        &
     &            +condensed_mix_ratio(l, i, ip_clcmp_cnv_ice)
            enddo
         else if (i_cloud_representation == ip_cloud_csiw) then
            do l=1, n_profile
               total_mix_ratio_st(l)                                    &
     &            =condensed_mix_ratio(l, i, ip_clcmp_st_water)
               total_mix_ratio_cnv(l)                                   &
     &            =condensed_mix_ratio(l, i, ip_clcmp_cnv_water)
            enddo
         endif
         do l=1, n_profile
            if (land_g(l)) then
               kparam=kparam_land
            else
               kparam=kparam_sea
            endif
            condensed_re(l, i, ip_clcmp_cnv_water)                      &
     &         =(3.0e+00*total_mix_ratio_cnv(l)*density_air(l, i)       &
     &         /(4.0e+00*pi*rho_water*kparam*n_drop(l, i)))             &
     &         **(1.0e+00/3.0e+00)
            condensed_re(l, i, ip_clcmp_st_water)                       &
     &         =(3.0e+00*total_mix_ratio_st(l)*density_air(l, i)        &
     &         /(4.0e+00*pi*rho_water*kparam*n_drop(l, i)))             &
     &         **(1.0e+00/3.0e+00)
         enddo
         do l=1, n_profile
            ntot_diag_g(l, i)=n_drop(l, i)*1.0e-06
            strat_lwc_diag_g(l, i)                                      &
     &         =total_mix_ratio_st(l)*density_air(l, i)*1.0e03
         enddo
      enddo
!
!     reset the effective radii for deep convective clouds.
      do i=n_layer+1-nclds, n_layer
         do l=1, n_profile
            if (land_g(l)) then
               if (cc_depth(l)  >   deep_convection_limit_land) then
                  condensed_re(l, i, ip_clcmp_cnv_water)=dconre_land
               else
                  if (condensed_re(l, i, ip_clcmp_cnv_water)            &
     &                                              >   dconre_land)    &
     &            condensed_re(l, i, ip_clcmp_cnv_water)=dconre_land
               endif
            else
               if (cc_depth(l)  >   deep_convection_limit_sea) then
                  condensed_re(l, i, ip_clcmp_cnv_water)=dconre_sea
               else
                  if (condensed_re(l, i, ip_clcmp_cnv_water)            &
     &                                              >   dconre_sea)     &
     &            condensed_re(l, i, ip_clcmp_cnv_water)=dconre_sea
               endif
            endif
         enddo
      enddo
!
!
!
      return
      END SUBROUTINE r2_re_mrf_umist
#endif
#endif
