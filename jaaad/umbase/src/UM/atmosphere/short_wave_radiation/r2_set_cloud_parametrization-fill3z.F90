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
      subroutine r2_set_cloud_parametrization(ierr, n_band              &
     &   , i_st_water, i_cnv_water, i_st_ice, i_cnv_ice                 &
     &   , l_drop_type, i_drop_parametrization                          &
     &   , n_drop_phf_term, drop_parameter_list                         &
     &   , drop_parm_min_dim, drop_parm_max_dim                         &
     &   , l_ice_type, i_ice_parametrization                            &
     &   , n_ice_phf_term, ice_parameter_list                           &
     &   , ice_parm_min_dim, ice_parm_max_dim                           &
     &   , i_condensed_param, condensed_n_phf, condensed_param_list     &
     &   , condensed_min_dim, condensed_max_dim                         &
     &   , nd_band, nd_drop_type, nd_ice_type, nd_cloud_parameter       &
     &   , nd_cloud_component                                           &
     &   )
!
!
!
      implicit none
!
!
!     comdecks included.
#include "cloud_component_pcf3z.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
!
!
!     dummy arguments:
!
      integer                                                           &
                !, intent(out)
     &     ierr
!             error flag
!
!     sizes of arrays:
      integer                                                           &
                !, intent(in)
     &     nd_band                                                      &
!             maximum number of spectral bands
     &   , nd_drop_type                                                 &
!             maximum number of types of droplets
     &   , nd_ice_type                                                  &
!             maximum number of types of ice crystals
     &   , nd_cloud_parameter                                           &
!             maximum number of parameters for clouds
     &   , nd_cloud_component
!             size allocated for components in clouds
!
      integer                                                           &
                !, intent(in)
     &     n_band
!             number of spectral bands
!
!     types of droplets and crystals:
      integer                                                           &
                !, intent(in)
     &     i_st_water                                                   &
!             type of water droplets in stratiform clouds
     &   , i_cnv_water                                                  &
!             type of water droplets in convective clouds
     &   , i_st_ice                                                     &
!             type of ice crystals in stratiform clouds
     &   , i_cnv_ice
!             type of ice crystals in convective clouds
!
      logical                                                           &
                !, intent(in)
     &     l_drop_type(nd_drop_type)                                    &
!             flags for types of droplet present
     &   , l_ice_type(nd_ice_type)
!             flags for types of ice crystal present
      integer                                                           &
                !, intent(in)
     &     i_drop_parametrization(nd_drop_type)                         &
!             parametrizations of types of droplets
     &   , n_drop_phf_term(nd_drop_type)                                &
!             Number of terms in the phase function for
!             droplets
     &   , i_ice_parametrization(nd_ice_type)                           &
!             parametrizations of types of ice crystals
     &   , n_ice_phf_term(nd_ice_type)
!             Number of terms in the phase function for
!             ice crystals
      real                                                              &
                !, intent(in)
     &     drop_parameter_list(nd_cloud_parameter                       &
     &        , nd_band, nd_drop_type)                                  &
!             parameters for optical parametrizations of droplets
     &   , drop_parm_min_dim(nd_drop_type)                              &
!             minimum size of droplets permitted in parametrizations
     &   , drop_parm_max_dim(nd_drop_type)                              &
!             maximum size of droplets permitted in parametrizations
     &   , ice_parameter_list(nd_cloud_parameter                        &
     &        , nd_band, nd_ice_type)                                   &
!             parameters for optical parametrizations of ice crystals
     &   , ice_parm_min_dim(nd_ice_type)                                &
!             minimum size of ice crystals permitted in parametrizations
     &   , ice_parm_max_dim(nd_ice_type)
!             maximum size of ice crystals permitted in parametrizations
!
      integer                                                           &
                !, intent(out)
     &     i_condensed_param(nd_cloud_component)                        &
!             types of parametrization used for condensed
!             components in clouds
     &   , condensed_n_phf(nd_cloud_component)
!             Number of terms in the phase function
      real                                                              &
                !, intent(out)
     &     condensed_param_list(nd_cloud_parameter                      &
     &        , nd_cloud_component, nd_band)                            &
!             coefficients for parametrization of condensed phases
     &   , condensed_min_dim(nd_cloud_component)                        &
!             minimum dimension of each condensed component
     &   , condensed_max_dim(nd_cloud_component)
!             maximum dimension of each condensed component
!
!
!     local variables:
      integer                                                           &
     &     i                                                            &
!             loop variable
     &   , j                                                            &
!             loop variable
     &   , i_scheme
!             parametrization scheme
!
!     functions called:
      integer                                                           &
     &     set_n_cloud_parameter
!             function to find number of parameters for clouds
      external                                                          &
     &     set_n_cloud_parameter
!
!
!
!     select parametrization for water in stratiform clouds:
!
      if ( (i_st_water <= nd_drop_type).and.                            &
     &     (l_drop_type(i_st_water)) ) then
         i_scheme=i_drop_parametrization(i_st_water)
         i_condensed_param(ip_clcmp_st_water)=i_scheme
         condensed_n_phf(ip_clcmp_st_water)=n_drop_phf_term(i_st_water)
         condensed_min_dim(ip_clcmp_st_water)                           &
     &      =drop_parm_min_dim(i_st_water)
         condensed_max_dim(ip_clcmp_st_water)                           &
     &      =drop_parm_max_dim(i_st_water)
      else
         write(iu_err, '(/a, /a)') '*** error: no data exist for type ' &
     &     , 'of droplet selected in stratiform water clouds.'
         ierr=i_err_fatal
         return
      endif
!
      do i=1, n_band
! DEPENDS ON: set_n_cloud_parameter
         do j=1, set_n_cloud_parameter(i_scheme                         &
     &      , ip_clcmp_st_water, condensed_n_phf(ip_clcmp_st_water))
            condensed_param_list(j, ip_clcmp_st_water, i)               &
     &         =drop_parameter_list(j, i, i_st_water)
         enddo
      enddo
!
!
!     select parametrization for water in convective clouds:
!
      if ( (i_cnv_water <= nd_drop_type).and.                           &
     &     (l_drop_type(i_cnv_water)) ) then
         i_scheme=i_drop_parametrization(i_cnv_water)
         i_condensed_param(ip_clcmp_cnv_water)=i_scheme
         condensed_n_phf(ip_clcmp_cnv_water)                            &
     &      =n_drop_phf_term(i_cnv_water)
         condensed_min_dim(ip_clcmp_cnv_water)                          &
     &      =drop_parm_min_dim(i_cnv_water)
         condensed_max_dim(ip_clcmp_cnv_water)                          &
     &      =drop_parm_max_dim(i_cnv_water)
      else
         write(iu_err, '(/a, /a)') '*** error: no data exist for type ' &
     &     , 'of crystal selected in convective water clouds.'
         ierr=i_err_fatal
         return
      endif
!
      do i=1, n_band
! DEPENDS ON: set_n_cloud_parameter
         do j=1, set_n_cloud_parameter(i_scheme                         &
     &      , ip_clcmp_cnv_water, condensed_n_phf(ip_clcmp_cnv_water))
            condensed_param_list(j, ip_clcmp_cnv_water, i)              &
     &         =drop_parameter_list(j, i, i_cnv_water)
         enddo
      enddo
!
!
!     select parametrization for ice in stratiform clouds:
!
      if ( (i_st_ice <= nd_ice_type).and.                               &
     &     (l_ice_type(i_st_ice)) ) then
         i_scheme=i_ice_parametrization(i_st_ice)
         i_condensed_param(ip_clcmp_st_ice)=i_scheme
         condensed_n_phf(ip_clcmp_st_ice)=n_ice_phf_term(i_st_ice)
         condensed_min_dim(ip_clcmp_st_ice)                             &
     &      =ice_parm_min_dim(i_st_ice)
         condensed_max_dim(ip_clcmp_st_ice)                             &
     &      =ice_parm_max_dim(i_st_ice)
      else
         write(iu_err, '(/a, /a)') '*** error: no data exist for type ' &
     &      , 'of crystal selected in stratiform ice clouds.'
         ierr=i_err_fatal
         return
      endif
!
      do i=1, n_band
! DEPENDS ON: set_n_cloud_parameter
         do j=1, set_n_cloud_parameter(i_scheme                         &
     &      , ip_clcmp_st_ice, condensed_n_phf(ip_clcmp_st_ice))
            condensed_param_list(j, ip_clcmp_st_ice, i)                 &
     &         =ice_parameter_list(j, i, i_st_ice)
         enddo
      enddo
!
!
!     select parametrization for ice in convective clouds:
!
      if ( (i_cnv_ice <= nd_ice_type).and.                              &
     &     (l_ice_type(i_cnv_ice)) ) then
         i_scheme=i_ice_parametrization(i_cnv_ice)
         i_condensed_param(ip_clcmp_cnv_ice)=i_scheme
         condensed_n_phf(ip_clcmp_cnv_ice)=n_ice_phf_term(i_cnv_ice)
         condensed_min_dim(ip_clcmp_cnv_ice)                            &
     &      =ice_parm_min_dim(i_cnv_ice)
         condensed_max_dim(ip_clcmp_cnv_ice)                            &
     &      =ice_parm_max_dim(i_cnv_ice)
      else
         write(iu_err, '(/a, /a)') '*** error: no data exist for type ' &
     &      , 'of crystal selected in convective ice clouds.'
         ierr=i_err_fatal
         return
      endif
!
      do i=1, n_band
! DEPENDS ON: set_n_cloud_parameter
         do j=1, set_n_cloud_parameter(i_scheme                         &
     &      , ip_clcmp_cnv_ice, condensed_n_phf(ip_clcmp_cnv_ice))
            condensed_param_list(j, ip_clcmp_cnv_ice, i)                &
     &         =ice_parameter_list(j, i, i_cnv_ice)
         enddo
      enddo
!
!
!
      return
      END SUBROUTINE r2_set_cloud_parametrization
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
