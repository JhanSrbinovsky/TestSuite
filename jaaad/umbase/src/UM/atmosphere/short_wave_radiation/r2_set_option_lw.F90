#if defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
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
      subroutine r2_set_option_lw(ierr                                  &
     &  , n_layer, n_band                                               &
     &  , isolir                                                        &
     &  , i_gas_overlap_in, i_gas_overlap                               &
     &  , l_aerosol_enabled, l_climat_aerosol                           &
     &  , l_use_sulpc_direct, l_use_soot_direct, l_use_biogenic         &
     &  , l_use_seasalt_direct, l_murk_rad, l_aerosol                   &
     &  , l_use_sulpc_indirect, l_aerosol_ccn, n_arcl_species           &
     &  , l_global_cloud_top, global_cloud_top, n_cloud_top_global      &
     &  , l_clear_olr, l_surf_down_clr, l_clear_hr, l_clear             &
     &  , i_angular_integration, i_solver, i_solver_clear               &
     &  , l_rescale, n_order_forward                                    &
     &  , i_truncation, ls_global_trunc, l_euler_trnf, euler_factor     &
     &  , i_scatter_method, i_scatter_method_band                       &
     &  , l_rad_tile                                                    &
     &  , weight_band                                                   &
     &  , nd_overlap_coeff, nd_2sg_profile, nd_layer_clr                &
     &  , nd_source_coeff, nd_max_order, nd_sph_coeff                   &
     &  , nd_region                                                     &
     &  , nd_profile, nd_band                                           &
     &  )
!
!
!
!     Modules included
      use tileid3z
!
      implicit none
!
!
!     Comdecks included
#include "solver_pcf3z.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
#include "spectral_region_pcf3z.h"
#include "angular_integration_pcf3z.h"
#include "sph_truncation_pcf3z.h"
!
!     Dummy variables:
!
!     Dimensions of arrays:
      integer                                                           &
                !, intent(in)
     &    nd_profile                                                    &
!           Maximum number of atmospheric profiles
     &  , nd_band
!           Maximum number of spectral bands
!
      integer                                                           &
                !, intent(out)
     &    ierr
!           Error flag
!
!     Actual sizes used:
      integer                                                           &
                !, intent(in)
     &    n_layer                                                       &
!           Number of atmospheric layers used in radiation
     &  , n_band
!           Number of spectral bands
!
!     Aerosol Flags:
      logical                                                           &
                !, intent(in)
     &    l_aerosol_enabled                                             &
!           Generic radiative switch for aerosol effects
     &  , l_climat_aerosol                                              &
!           Flag for climatological aerosols
     &  , l_use_sulpc_direct                                            &
!           Flag to include the direct effect of sulphate aerosols
     &  , l_use_soot_direct                                             &
!           Flag to include the direct effect of soot aerosols
     &  , l_use_sulpc_indirect                                          &
!           Flag to include the indirect effect of sulphate aerosols
     &  , l_use_seasalt_direct                                          &
!           Flag to include the direct effect of seasalt aerosols
     &  , l_use_biogenic                                                &
!           Flag to include the direct effect of biogenic aerosols
     &  , l_murk_rad
!           Flag to include urban aerosols

      integer n_arcl_species
!           Number of species from the NWP aerosol climatology.
!           Zero is the NWP climatology is not used.

      logical                                                           &
                !, intent(out)
     &    l_aerosol                                                     &
!           Flag for direct aerosol effects passed to the radiation
!           code
     &  , l_aerosol_ccn
!           Flag for indirect aerosol effects passed to the radiation
!           code
!
      integer                                                           &
                 !, intent(in)
     &    i_angular_integration                                         &
!           Method of angular integration
     &  , i_solver                                                      &
!           Solver selected
     &  , global_cloud_top                                              &
!           Cloud top in ascending order
     &  , i_gas_overlap_in                                              &
!           Overall option for gaseous overlaps
     &  , isolir                                                        &
!           Spectral region (technically redundant, but I'd eventually
!           like to unify with the SW)
     &  , i_scatter_method
!           Treatment of scattering supplied in the controlling list
      logical                                                           &
                 !, intent(in)
     &    l_rescale                                                     &
!           Rescaling flag
     &  , l_euler_trnf
!           Flag to apply Euler's transformation to alternating series
      logical                                                           &
                !, intent(in)
     &    l_global_cloud_top                                            &
!           Flag to use a global cloud-top
     &  , l_clear_olr                                                   &
!           Flag for clear-sky outgoing longwave flux
     &  , l_surf_down_clr                                               &
!           Flag for clear-sky downward surface flux
     &  , l_clear_hr
!           Flag for clear-sky heating rates
!
!
!
!     Dimensions defined by the routine
      integer                                                           &
                !, intent(out)
     &    nd_source_coeff                                               &
!           Size allocated for two-stream source coefficients
     &  , nd_max_order                                                  &
!           Size allocated for orders of spherical harmonics
     &  , nd_sph_coeff                                                  &
!           Size allocated for coefficients of spherical harmonics
!           (includes polar and azimuthal orders)
     &  , nd_overlap_coeff                                              &
!           Size allocated for cloud overlap coefficients
     &  , nd_region                                                     &
!           Size allocated for aggregated cloudy regions
     &  , nd_2sg_profile                                                &
!           Size allocated for profiles in two-stream flux arrays
     &  , nd_layer_clr
!           Size allocated for totally clear layers
      integer                                                           &
                 !, intent(out)
     &    i_solver_clear                                                &
!           Clear-sky solver
     &  , n_cloud_top_global                                            &
!           Inverted global value for the cloud top
     &  , i_gas_overlap(n_band)                                         &
!           Gaseous overlaps in each band
     &  , ls_global_trunc                                               &
!           Order of global truncation
     &  , i_truncation                                                  &
!           Type of spherical truncation
     &  , n_order_forward
!           Order of term in phase function used to define forward
!           scattering fraction
      integer                                                           &
                !, intent(out)
     &    i_scatter_method_band(nd_band)
!           Treatment of scattering in each band
      real                                                              &
                !, intent(out)
     &    euler_factor                                                  &
!           Weighting applied to the last term of alternating series
     &  , weight_band(nd_band)
!           Weightings applied to each band
      logical                                                           &
                !, intent(out)
     &    l_clear                                                       &
!           Flag for clear-sky calculations
     &  , l_rad_tile
!           Flag to enable surface tiling
!
!
!
!     Local variables.
      integer                                                           &
     &    i
!           Loop variable
!
!
!
!     Decide on the final optins for aerosols:
      l_aerosol=l_aerosol_enabled.and.                                  &
     &  (l_use_sulpc_direct.or.                                         &
     &   l_use_soot_direct.or.                                          &
     &   l_use_seasalt_direct.or.                                       &
     &   l_use_biogenic.or.                                             &
     &   l_murk_rad.or.                                                 &
     &   l_climat_aerosol.or.                                           &
     &   (n_arcl_species > 0))
!
!     Whilst l_aerosol_ccn is a generic flag for determining CCN
!     from aerosol, the view is currently taken that sulphate aerosols
!     must be included with all indirect effects, other aerosols
!     being additional, so l_aerosol_ccn is assigned solely from
!     l_use_sulpc_indirect.
      l_aerosol_ccn=l_use_sulpc_indirect
!
      if (i_angular_integration == ip_two_stream) then
!
        if (isolir == ip_solar) then
          if ( (i_solver /= ip_solver_pentadiagonal).and.               &
     &         (i_solver /= ip_solver_mix_11).and.                      &
     &         (i_solver /= ip_solver_mix_direct).and.                  &
     &         (i_solver /= ip_solver_mix_direct_hogan).and.            &
     &         (i_solver /= ip_solver_homogen_direct).and.              &
     &         (i_solver /= ip_solver_triple).and.                      &
     &         (i_solver /= ip_solver_triple_hogan)                     &
     &      ) then
             write(iu_err, '(/a, /a)')                                  &
     &          '*** error: an invalid two-stream solver has been '     &
     &        , 'selected in the shortwave region.'
             ierr=i_err_fatal
             return
          endif
!
        else if (isolir == ip_infra_red) then
!
          if ( (i_solver /= ip_solver_pentadiagonal).and.               &
     &         (i_solver /= ip_solver_mix_11).and.                      &
     &         (i_solver /= ip_solver_mix_app_scat).and.                &
     &         (i_solver /= ip_solver_mix_direct).and.                  &
     &         (i_solver /= ip_solver_mix_direct_hogan).and.            &
     &         (i_solver /= ip_solver_homogen_direct).and.              &
     &         (i_solver /= ip_solver_triple).and.                      &
     &         (i_solver /= ip_solver_triple_hogan).and.                &
     &         (i_solver /= ip_solver_triple_app_scat)                  &
     &      ) then
            write(iu_err, '(/a, /a)')                                   &
     &        '*** error: an invalid solver has been selected '         &
     &        , 'in the longwave region.'
            ierr=i_err_fatal
            return
          endif
!
        endif
!
        nd_2sg_profile=nd_profile
        nd_source_coeff=2
        nd_max_order=1
        nd_sph_coeff=0
        if ( (i_solver == ip_solver_triple).or.                         &
     &       (i_solver == ip_solver_triple_hogan).or.                   &
     &       (i_solver == ip_solver_triple_app_scat) ) then
          nd_overlap_coeff=18
        else if ( (i_solver == ip_solver_mix_direct).or.                &
     &       (i_solver == ip_solver_mix_direct_hogan).or.               &
     &       (i_solver == ip_solver_mix_11) ) then
          nd_overlap_coeff=8
        else
          nd_overlap_coeff=0
        endif
!
        if ( (i_solver == ip_solver_triple).or.                         &
     &       (i_solver == ip_solver_triple_hogan).or.                   &
     &       (i_solver == ip_solver_triple_app_scat) ) then
          nd_region=3
        else
          nd_region=2
        endif
!
        if (l_rescale) n_order_forward=2
!
!       Set clear-sky calculations.
        l_clear=l_clear_olr.or.                                         &
     &          l_surf_down_clr.or.                                     &
     &          l_clear_hr
!
        if (l_clear) then
!
!         Select a clear-sky solver to match the main solver.
          if (i_solver == ip_solver_pentadiagonal) then
            i_solver_clear=ip_solver_pentadiagonal
          else if (i_solver == ip_solver_mix_11) then
            i_solver_clear=ip_solver_pentadiagonal
          else if (i_solver == ip_solver_mix_direct) then
            i_solver_clear=ip_solver_homogen_direct
          else if (i_solver == ip_solver_mix_direct_hogan) then
            i_solver_clear=ip_solver_homogen_direct
          else if (i_solver == ip_solver_homogen_direct) then
            i_solver_clear=ip_solver_homogen_direct
          else if (i_solver == ip_solver_triple) then
            i_solver_clear=ip_solver_homogen_direct
          else if (i_solver == ip_solver_triple_hogan) then
            i_solver_clear=ip_solver_homogen_direct
          else if (i_solver == ip_solver_triple_app_scat) then
            i_solver_clear=ip_solver_homogen_direct
          else
            write(iu_err, '(/a, /a)')                                   &
     &        '*** error: no clear-sky counterpart has been '           &
     &        , 'specified for this lw solver.'
            ierr=i_err_fatal
            return
          endif
!
        endif
!
!       We permit tiling of sea-ice points only with the two-stream
!       option at present.
        l_rad_tile=.true.
!
      else if (i_angular_integration == ip_spherical_harmonic) then
!
        nd_2sg_profile=1
        nd_source_coeff=0
        nd_overlap_coeff=0
        nd_region=0
        nd_max_order=ls_global_trunc+2
        if (i_truncation == ip_trunc_triangular) then
          nd_sph_coeff=(ls_global_trunc+3)*(ls_global_trunc+4)/2
        else if (i_truncation == ip_trunc_azim_sym) then
          nd_sph_coeff=ls_global_trunc+2
        else
          write(iu_err, '(/a)')                                         &
     &      'error: illegal truncation'
          ierr=i_err_fatal
          return
        endif
!
        if (l_rescale) n_order_forward=ls_global_trunc+1
!
!       As currently implemented, Euler's transformation is applied
!       only in its most basic form, adding just half of the last
!       term in an alternating series.
        if (l_euler_trnf) then
          euler_factor=0.5
        else
          euler_factor=1.0
        endif
!
!       Clear-sky fluxes are not available from the spherical harmonic
!       code in the same call as cloudy fluxes yet. If required, they
!       should be diagnosed by using a separate call to the code with
!       clouds switched off.
        l_clear=.false.
        if (l_clear_olr     .OR.                                        &
     &      l_surf_down_clr .OR.                                        &
     &      l_clear_hr) then
          ierr=i_err_fatal
! DEPENDS ON: ereport
          call Ereport('set_option_lw', ierr                            &
     &      , 'Clear-sky fluxes not directly available in harmonics')
        endif
!
!       We permit tiling of sea-ice points only with the two-stream
!       option at present.
        l_rad_tile=.false.
!
      endif
!
!
!
!     Set properties for individual bands.
      do i=1, n_band
        weight_band(i)=1.0
        i_gas_overlap(i)=i_gas_overlap_in
!       Extend the treatment of scattering from the control structure
!       to each band.
        i_scatter_method_band(i)=i_scatter_method
      enddo
!
!
!     Invert the topmost cloudy layer if using a global value.
      if (l_global_cloud_top) then
        n_cloud_top_global=n_layer+1-global_cloud_top
        nd_layer_clr=n_cloud_top_global-1
      endif
!
!
!
      return
      END SUBROUTINE r2_set_option_lw
#endif
