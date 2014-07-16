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
! Description of code:
!   FORTRAN 90
!
!- ---------------------------------------------------------------------
      subroutine r2_lwrad3z(ierr                                        &
!                       Gaseous mixing ratios
     &  , h2o, co2, o3                                                  &
     &  , co2_dim1, co2_dim2, co2_3d, l_co2_3d                          &
     &  , n2o_mix_ratio, ch4_mix_ratio                                  &
     &  , cfc11_mix_ratio, cfc12_mix_ratio, cfc113_mix_ratio            &
     &  , hcfc22_mix_ratio, hfc125_mix_ratio, hfc134a_mix_ratio         &
!                       Thermodynamic variables
     &  , tac,  tstar,  tstar_solid, tstar_sea, l_ctile, pstar          &
     &  , p_layer_boundaries                                            &
     &  , p_layer_centres                                               &
     &  , height_theta                                                  &
     &  , height_rho                                                    &
!                       Options for treating clouds
     &  , global_cloud_top, l_inhom_cloud, inhom_cloud                  &
     &  , dp_corr_strat, dp_corr_conv                                   &
!                       Stratiform cloud fields
     &  , l_cloud_water_partition                                       &
     &  , lca_area, lca_bulk, lccwc1, lccwc2                            &
!                       Convective cloud fields
     &  , cca, cccwp, ccb, cct                                          &
!                       Surface fields
     &  , land, flandg, ice_fraction                                    &
     &  , lying_snow                                                    &
!                       Aerosol fields
     &  , l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero        &
     &  , bl_depth, n_levels_bl                                         &
     &  , l_use_sulpc_direct, l_use_sulpc_indirect                      &
     &  , sulp_dim1,sulp_dim2                                           &
     &  , accum_sulphate, aitken_sulphate, diss_sulphate                &
     &  , sea_salt_film, sea_salt_jet, l_use_seasalt_indirect           &
     &  , l_use_seasalt_direct, salt_dim_a, salt_dim_b                  &
     &  , l_use_soot_direct, soot_dim1, soot_dim2                       &
     &  , fresh_soot, aged_soot                                         &
     &  , l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic        &
     &  , l_use_arcl, arcl_dim1, arcl_dim2, n_arcl_species              &
     &  , n_arcl_compnts, i_arcl_compnts, arcl                          &
     &  , aero_meso, l_murk_rad, Ntot_land, Ntot_sea                    &
!                       Level of tropopause
     &  , trindx                                                        &
!                       Spectrum
     &  , lw_spectrum                                                   &
!                       Algorithmic options
     &  , lw_control                                                    &
     &  , pts, l_mod_k_flux, l_scale_inc, list, i_call                  &
#if defined(RAD_DBG)
     &  , i_segment                                                     &
#endif
!                       Satellite viewing geometry
     &  , n_viewing_direction                                           &
     &  , viewing_direction1, viewing_direction2                        &
     &  , n_viewing_level, viewing_level                                &
!                       Diagnostics
     &  , n_channel, map_channel                                        &
     &  , LW_diag, row_list, col_list                                   &
!                       Physical dimensions
     &  , n_points, nlevs, n_layer, nclds                               &
     &  , nwet, nozone, row_length, rows, nd_field                      &
     &  , nd_field_flux_diag, nd_field_rad_diag                         &
     &  , nd_profile, nd_layer, nd_column, n_cca_lev, nd_channel        &
     &  , nd_flux_profile, nd_radiance_profile                          &
     &  , nd_viewing_level, nd_direction                                &
     &  , nd_cloud_component, nd_cloud_type                             &
     &  , nd_brdf_basis_fnc, nd_brdf_trunc                              &
     &  , nd_point_tile, nd_tile, id_ct                                 &
!                       Output fields
     &  , olr, top_absorption, lwsea, lwout                             &
           ! Variables needed to calculate layer masses
     &   , rho_r2, r_rho_levels, r_theta_levels                         &
     &   , q, qcl, qcf, qcf2, qrain, qgraup                             &
     &   , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
     &  )
!
!
!
!     Modules included.
      use dec_spec
      use control_struc
      use tileid3z

      Use lwrdiag_mod: Only                                             &
          StrLWDiag
!
      implicit none
!
!
!     Comdecks included
#include "c_0_dg_c.h"
#include "c_kinds.h"
!
!     Internal dimensions of the code
!     Spectral regions
#include "spectral_region_pcf3z.h"
!     Methods of integration
#include "angular_integration_pcf3z.h"
!     Mode of operation
#include "sph_mode_pcf3z.h"
!     Algebraic solvers
#include "solver_pcf3z.h"
!     Methods of scattering
#include "scatter_method_pcf3z.h"
!     Physical constants
#include "physical_constants_0_ccf3z.h"
!     Mathematical constants
#include "c_pi.h"
!     Unit numbers for printed output
#include "def_std_io_icf3z.h"
!     Error flags
#include "error_pcf3z.h"
!    Components of clouds for two-stream radiation code.
#include "cloud_component_pcf3z.h"
!    Cloud types for two-stream radiation code.
#include "cloud_type_pcf3z.h"
!
!
!     Dummy arguments
!
      integer                                                           &
                !, intent(out)
     &    ierr
!           Error flag
!
!     Dimensions of arrays:
      integer, intent(in) :: row_length
!                              length of rows on each domain
      integer, intent(in) :: rows
!                              number of rows in the domain
      integer                                                           &
                !, intent(in)
     &    nd_field                                                      &
!           Field size in calling program
     &  , nd_field_flux_diag                                            &
!           Field size for flux diagnostics
     &  , nd_field_rad_diag                                             &
!           Field size for radiance diagnostics
     &  , nd_profile                                                    &
!           Size of array of grid-points
     &  , nd_layer                                                      &
!           Array sizes for layers
     &  , nd_column
!           Number of columns per point
      integer                                                           &
                !, intent(in)
     &    nd_channel                                                    &
!           Size allocated for spectral channels
     &  , nd_flux_profile                                               &
!           Size allocated for grid-points where fluxes are
!           diagnosed
     &  , nd_radiance_profile                                           &
!           Number of points where radiances are required
     &  , nd_viewing_level                                              &
!           Number of levels where radiances are required
     &  , nd_direction                                                  &
!           Number of directions in which radiances are calculated
     &  , nd_cloud_component                                            &
!           Number of components permitted in clouds
     &  , nd_cloud_type                                                 &
!           Number of permitted types of cloud
     &  , nd_brdf_basis_fnc                                             &
!           Maximum permitted number of BRDF basis functions
     &  , nd_brdf_trunc                                                 &
!           Maximum permitted order of truncation for BRDFs
     &  , nd_point_tile                                                 &
!           Size allocated for points where the surface is tiled
     &  , nd_tile
!           Size allocated for surface tiles
!
!     Actual sizes used:
      integer                                                           &
                !, intent(in)
     &    n_points                                                      &
!           Number of points
     &  , nwet                                                          &
!           Number of wet levels
     &  , nozone                                                        &
!           Number of levels with ozone
     &  , nlevs                                                         &
!           Number of layers in the main model
     &  , n_layer                                                       &
!           Number of layers seen in the radiation scheme
     &  , nclds                                                         &
!           Number of cloudy levels
     &  , n_levels_bl                                                   &
!           Number of layers occupied by boundary-layer aerosol
!           if l_clim_aero_hgt is false.
     &  , n_cca_lev
!
!     Spectral data:
      type (spectrum) lw_spectrum
!
!     Control data:
      type (control_option) lw_control
!
      logical                                                           &
     &    l_scale_inc
!           Flag for scaling of heating rates to increments
      INTEGER :: i_call
!        Number of calls
!
#if defined(RAD_DBG)
      INTEGER :: i_segment
!       Segment in call
#endif
!
!
!
!     Gaseous mixing ratios:
      real                                                              &
                !, intent(in)
     &    h2o(nd_field, nwet)                                           &
!           Mass mixing ratio of water
     &  , co2                                                           &
!           Mass mixing ratio of CO2
     &  , o3(nd_field, nozone)                                          &
!           Mass mixing ratios of ozone
     &  , n2o_mix_ratio                                                 &
!           Mass mixing ratio of nitrous oxide
     &  , ch4_mix_ratio                                                 &
!           Mass mixing ratio of methane
     &  , cfc11_mix_ratio                                               &
!           Mass mixing ratio of CFC11
     &  , cfc12_mix_ratio                                               &
!           Mass mixing ratio of CFC12
     &  , cfc113_mix_ratio                                              &
!           Mass mixing ratio of CFC113
     &  , hcfc22_mix_ratio                                              &
!           Mass mixing ratio of HCFC22
     &  , hfc125_mix_ratio                                              &
!           Mass mixing ratio of HFC125
     &  , hfc134a_mix_ratio
!           Mass mixing ratio of HFC134a
!
!     General atmospheric properties:
      real                                                              &
                !, intent(in)
     &    p_layer_boundaries(nd_field,0:nlevs)                          &
!             pressure at boundaries of layers
     &  , p_layer_centres(nd_field,0:nlevs)                             &
!             pressure at centres of layers
     &  , height_theta(nd_field, 0:nlevs)                               &
     &  , height_rho(nd_field, nlevs)                                   &
     &  , tac(nd_field, nlevs)
!           Temperatures at centres of layers

      Logical, intent(in)::                                             &
     &     l_mcr_qcf2                                                   &
                          ! Use second ice category
     &,    l_mcr_qrain                                                  &
                          ! Use prognostic rain
     &,    l_mcr_qgraup                                                 &
                          ! Use graupel
     &,    l_mixing_ratio ! Use mixing ratios in layer mass calculation
!
!     Options for treating clouds
      integer                                                           &
                !, intent(in)
     &    global_cloud_top
!           Global topmost cloudy layer
      logical                                                           &
                !, intent(in)
     &    l_inhom_cloud
!           Flag to use scaling factors for inhomogeneous cloud
      real                                                              &
                !, intent(in)
     &    inhom_cloud(nd_cloud_component)                               &
!           Scaling factors for inhomogeneous cloud
     &  , dp_corr_strat                                                 &
!           Decorrelation pressure scale for large scale cloud
     &  , dp_corr_conv
!           Decorrelation pressure scale for convective cloud
!
!     Properties of stratiform clouds:
      logical                                                           &
                !, intent(in)
     &    l_cloud_water_partition
!           Flag to use prognostic cloud ice contents
      real                                                              &
                !, intent(in)
     &    lccwc1(nd_field, nclds+1/(nclds+1))                           &
!           Liquid water contents (these are not used directly in
!           the radiation: the total condensed water content is
!           repartitioned using focwwil).
     &  , lccwc2(nd_field, nclds+1/(nclds+1))                           &
!           Ice water contents (these are not used directly in
!           The radiation: the total condensed water content is
!           Repartitioned using focwwil).
     &  , lca_area(nd_field, nclds+1/(nclds+1))                         &
!           Area fractions of layer clouds outside convective towers
     &  , lca_bulk(nd_field, nclds+1/(nclds+1))
!           Bulk fractions of layer clouds outside convective towers
!
!     Properties of convective clouds:
      integer                                                           &
                !, intent(in)
     &    ccb(nd_field)                                                 &
!           Base of convective cloud
     &  , cct(nd_field)
!           Top of convective cloud
      real                                                              &
                !, intent(in)
     &    cccwp(nd_field)                                               &
!           Water path of convective cloud
     &  , cca(nd_field,n_cca_lev)
!           Fraction of grid-box covered by convective cloud

!     Aerosols:
      logical                                                           &
                !, intent(in)
     &    l_climat_aerosol                                              &
!           Flag for climatological aerosol
     &  , l_clim_aero_hgt                                               &
!           Flag to use the depth of the boundary layer to set
!           the climatological aerosol
     &   , L_HadGEM1_Clim_Aero                                          &
!           Flag to use HadGEM1 setting for climatological aerosols
     &  , l_murk_rad
!           Flag for mesoscale model aerosol
      logical                                                           &
                !, intent(in)
     &    l_use_sulpc_direct                                            &
!           Flag to use sulphur cycle for direct effect
     &  , l_use_sulpc_indirect                                          &
!           Flag to use sulphur cycle for indirect effect
     &  , l_use_soot_direct                                             &
!           Use direct rad. effect of soot aerosol
     &  , l_use_seasalt_indirect                                        &
!           Flag to use sea-salt for indirect effect
     &  , l_use_seasalt_direct                                          &
!           Flag to use sea-salt for direct effect
     &  , l_use_biogenic
!           Flag to use biogenic for direct effect
      integer                                                           &
                !,intent (in)
     &    sulp_dim1,sulp_dim2                                           &
!           Dimensions for _sulphate arrays, (P_FIELD,P_LEVELS or 1,1)
     &  , soot_dim1, soot_dim2                                          &
!           dimensions for soot arrays (P_FIELD,P_LEVELS or 1,1)
     &  , salt_dim_a, salt_dim_b                                        &
!           dimensions for salt arrays on input (salt_dim_a=p_field
!           and salt_dim_b=p_levels, or else 1,1)
     &  , salt_dim_ind_a, salt_dim_ind_b                                &
!           dimensions for sea-salt arrays passed down to
!           r2_set_cloud_field if indirect effect required.
     &  , salt_dim_dir_a, salt_dim_dir_b                                &
!           dimensions for sea-salt arrays passed down to
!           r2_set_aerosol_field if direct effect required.
     &  , biogenic_dim1, biogenic_dim2
!           dimensions for biogenic array passed down to
!           r2_set_aerosol_field if direct effect required.
      real                                                              &
                !, intent(in)
     &    accum_sulphate(sulp_dim1, sulp_dim2)                          &
!           Mass mixing ratio of accumulation mode aerosol
     &  , aitken_sulphate(sulp_dim1, sulp_dim2)                         &
!           Mass mixing ratio of aitken mode aerosol
     &  , diss_sulphate(sulp_dim1, sulp_dim2)                           &
!           Mixing ratio of dissolved sulphate
     &  , fresh_soot(soot_dim1, soot_dim2)                              &
!           Mixing ratios of fresh soot
     &  , aged_soot(soot_dim1, soot_dim2)                               &
!           Mixing ratios of aged soot
     &  , sea_salt_film(salt_dim_a, salt_dim_b)                         &
!           Number concentration of film-mode sea-salt aerosol
     &  , sea_salt_jet(salt_dim_a, salt_dim_b)                          &
!           Number concentration of jet-mode sea-salt aerosol
     &  , biogenic(biogenic_dim1, biogenic_dim2)                        &
!           Mixing ratios of biogenic aerosol
     &  , aero_meso(nd_field, nlevs)
!           Mixing ratio of 'urban' aerosol of mesoscale model
!
!     Aerosol climatology for NWP
#include "arcl_dim.h"

      ! Number of requested species within the climatology
      integer n_arcl_species
      
      ! Corresponding number of requested components
      integer n_arcl_compnts
      
      ! Model switch for each species
      logical l_use_arcl(NPD_ARCL_SPECIES)
      
      ! Index of each component
      integer i_arcl_compnts(NPD_ARCL_COMPNTS)
      
      ! Array dimensions
      integer                                                           &
     &        arcl_dim1                                                 &
     &   ,    arcl_dim2
     
      ! Mass-mixing ratios 
      real                                                              &
     &        arcl(arcl_dim1, arcl_dim2, n_arcl_compnts)
!
!     Carbon cycle:
      logical   l_co2_3d    !controls use of 3d co2 field
      integer                                                           &
                !, intent(in)
     &    co2_dim1, co2_dim2
!           Dimensions for CO2 array, (P_FIELD,P_LEVELS or 1,1)
      real                                                              &
                !, intent(in)
     &    co2_3d(co2_dim1, co2_dim2)
!           Mass mixing ratio of carbon dioxide
!     Surface fields:
      logical                                                           &
                !, intent(in)
     &    land(nd_field)                                                &
!           Land mask (true if land fraction >0.5)
     &  , l_ctile
!           coastal tiling switch
      real                                                              &
                !, intent(in)
     &    flandg(nd_field)
!           land fraction in grid box
      real                                                              &
                !, intent(in)
     &    pstar(nd_field)                                               &
!           Surface pressures
     &  , tstar(nd_field)                                               &
!           Surface temperatures
     &  , tstar_solid(nd_field)                                         &
!           solid surface temperature
!           (i.e. areal mean of land and sea-ice)
     &  , tstar_sea(nd_field)                                           &
!           open sea surface temperatures
     &  , ice_fraction(nd_field)                                        &
!           Sea ice fraction of sea portion of grid box
     &  , lying_snow(nd_field)                                          &
!           Mass loading of lying snow
     &  , bl_depth(nd_field)
!           depth of the boundary layer
!
      real                                                              &
                !, intent(in)
     &     Ntot_land                                                    &
!           number of droplets over land / m-3
     &   , Ntot_sea
!           number of droplets over sea / m-3
!
!                       Level of tropopause
      integer, intent(in) ::                                            &
     &    trindx(nd_field)
!           The layer boundary of the tropopause
!
!     increment of time:
      real                                                              &
                !, intent(in)
     &    pts
!           Time increment
!
!     Use modulus of fluxes to remove negative effective extinctions
      logical, intent(in) :: l_mod_k_flux

      integer, intent(in) :: list(nd_field)
!                              list of points where radiation is to be
!                              calculated
!
!     Satellite viewing geometry
      INTEGER, Intent(IN) :: n_channel
!           Number of channels calculated simultaneously in one
!           call to the radiation code
      INTEGER, Intent(IN) :: map_channel(lw_spectrum%npd_band)
!           Mapping of bands in the spectral file to channels in the
!           diagnostic output
      INTEGER, Intent(IN) :: n_viewing_direction
!           Number of viewing directions
      INTEGER, Intent(IN) :: n_viewing_level
!           Number of levels where the radiance is calculated
      REAL, Intent(IN) :: viewing_direction1(nd_field, nd_direction)
!           Satellite viewing direction (Zenith Angle)
      REAL, Intent(IN) :: viewing_direction2(nd_field, nd_direction)
!           Satellite viewing directions (Azimuthal Angle)
      REAL, Intent(IN) :: viewing_level(nd_viewing_level)
!           Levels where radiances are calculated
!
!
!     Information for the calculation of layer masses
      Real, intent(in)::                                                &
     &  rho_r2(nd_field,nlevs)                                          &
                                ! Air density*radius of earth**2 / kg m-1
     &, r_rho_levels(nd_field,nlevs)                                    &
                                      ! Height of rho levels / m
     &, r_theta_levels(nd_field,0:nlevs)                                &
                                           ! Height of theta levels / m
     &, q(nd_field,nwet)                                                &
                                ! Water vapour mixing ratio / kg kg-1
     &, qcl(nd_field,nwet)                                              &
                                ! Liquid water mixing ratio / kg kg-1
     &, qcf(nd_field,nwet)                                              &
                                ! Ice mixing ratio / kg kg-1
     &, qcf2(nd_field,nwet)                                             &
                                ! Second ice category mr / kg kg-1
     &, qrain(nd_field,nwet)                                            &
                                ! Rain mixing ratio / kg kg-1
     &, qgraup(nd_field,nwet)  ! Graupel mixing ratio / kg kg-1

!     Calculated fluxes:
      real                                                              &
                !, intent(out)
     &    olr(nd_field)                                                 &
!           Net outgoing radiation
     &  , top_absorption(nd_field)                                      &
!           Absorption in the extra radiative layer at the top
!           of the model
     &  , lwout(nd_field, nlevs+1)                                      &
!           Net downward fluxes or heating rates
     &  , lwsea(nd_field)
!           Sea-surface components of flux
!
!
!
!     Diagnostics:
      type (strlwdiag) :: lw_diag
!
      integer, intent(in) :: row_list(nd_field)
!                              list of row indices of points
!                              to be treated
      integer, intent(in) :: col_list(nd_field)
!                              list of column indices of points
!                              to be treated
!
!
!
!     Local variables.
!
!     Locally allocated dimensions
      integer                                                           &
     &    nd_source_coeff                                               &
!           Size allocated for source coefficients
     &  , nd_layer_clr                                                  &
!           Size allocated for totally clear layers
     &  , nd_2sg_profile                                                &
!           Size allocated for profiles in two-stream or Gaussian
!           runs
     &  , nd_region                                                     &
!           Size allocated for cloudy regions
     &  , nd_overlap_coeff                                              &
!           Size allocated for cloudy overlap coefficients
     &  , nd_max_order                                                  &
!           Size allocated for orders of spherical harmonics
     &  , nd_sph_coeff                                                  &
!           Size allocated for coefficients of spherical harmonics
     &  , id_ct
!           Top level in arrays of cloud properties
!
      integer                                                           &
     &    i                                                             &
!           Loop variable
     &  , j                                                             &
!           Loop variable
     &  , k                                                             &
!           Loop variable
     &  , l                                                             &
!           Loop variable
     &  , ll                                                            &
!           Loop variable
     &  , lll                                                           &
!           Loop variable
     &  , ic
!           Loop variable
      logical                                                           &
     &    l_clear
!           Calculate clear-sky fields
!     Flags for processes actually enabled.
      integer                                                           &
     &    i_solver_clear                                                &
!           Solver for clear-sky fluxes
     &  , i_gas_overlap(lw_spectrum%npd_band)                           &
!           Overlaps in each band
     &  , n_order_forward
!           Order of term used to define the
!           forward scattering fraction
!
!     General atmospheric properties:
      real                                                              &
     &    d_mass(nd_profile, nd_layer)                                  &
!           Mass thicknesses of layers
     &  , p(nd_profile, nd_layer)                                       &
!           Pressure field
     &  , t(nd_profile, nd_layer)                                       &
!           Temperature field
     &  , t_bdy(nd_profile, 0: nd_layer)                                &
!           Temperature field at boundaries
     &  , gas_mix_ratio(nd_profile, nd_layer, lw_spectrum%npd_species)
!           Mass fractions of gases
!
      real :: layer_heat_capacity(nd_profile, nd_layer)
!           Specific heat capacity of layer * d_mass

!     Surface fields:
      Logical :: Land_g(nd_profile)
!           Gathered land-surface mask
      Real :: ice_fraction_g(nd_profile)
!           Gathered ice fraction
      integer                                                           &
     &    n_brdf_basis_fnc
!           Number of BRDF basis functions
      real                                                              &
     &    rho_alb(nd_profile, nd_brdf_basis_fnc, lw_spectrum%npd_band)  &
!           Weights for BRDF basis functions
     &  , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                  &
     &      , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)                     &
!           Array of BRDF basis terms
     &  , t_surface(nd_profile)                                         &
!           Gathered temperature of surface
     &  , t_solid(nd_profile)                                           &
!             Gathered temperature of solid surface
     &  , t_sea(nd_profile)                                             &
!             GATHERED OPEN SEA TEMPERATURE
     &  , flandg_g(nd_profile)
!           Gathered land fraction
!
!     Arrays related to tiling of the surface
      logical                                                           &
     &    l_rad_tile
!           Local to allow tiling within the radiation scheme: this is
!           used to handle all surface heterogeneities, whether or not
!           explicit coastal tiling is required (i.e. even when a grid
!           box may be only land or sea, this facility is still used to
!           represent the split between open sea and sea-ice.)
      integer                                                           &
     &    n_point_tile                                                  &
!           Number of points to tile
     &  , n_tile                                                        &
!           Number of tiles used
     &  , list_tile(nd_point_tile)                                      &
!           List of points with surface tiling indexed over the list
!           of points where radiances are calculated
     &  , list_tile_outer(nd_point_tile)                                &
!           List of points with surface tiling indexed over the full
!           list of points where this routine has been called
     &  , index_tile(npd_tile_type)
!           The indexing number of tiles of the given type
      real                                                              &
     &    rho_alb_tile(nd_point_tile, nd_brdf_basis_fnc                 &
     &      , nd_tile, lw_spectrum%npd_band)                            &
!           Weights for the basis functions of the BRDFs
!           at the tiled points
     &  , t_tile(nd_point_tile, nd_tile)                                &
!           Local surface temperatures on individual tiles
     &  , frac_tile(nd_point_tile, nd_tile)
!           Fractions of ech tiled grid-point occupied by tiles of
!           the relevant types
!
!     Cloudy properties:
      integer                                                           &
     &    n_condensed                                                   &
!           Number of condensed phases
     &  , type_condensed(nd_cloud_component)                            &
!           Types of condensed components
     &  , i_condensed_param(nd_cloud_component)                         &
!           Parametrization schemes for components
     &  , condensed_n_phf(nd_cloud_component)                           &
!           Number of terms in the phase functions
     &  , n_cloud_top_global                                            &
!           Inverted global topmost cloudy layer
     &  , n_cloud_type
!           Number of types of clouds
      real                                                              &
     &    condensed_param_list(lw_spectrum%npd_cloud_parameter          &
     &      , nd_cloud_component, lw_spectrum%npd_band)                 &
!           Parameters for condensed phases
     &  , condensed_dim_char(nd_profile, id_ct: nd_layer                &
     &      , nd_cloud_component)                                       &
!           Characteristic dimensions of condensed species
     &  , condensed_mix_ratio(nd_profile, id_ct: nd_layer               &
     &      , nd_cloud_component)                                       &
!           Mass fractions of liquid water
     &  , w_cloud(nd_profile, id_ct: nd_layer)                          &
!           Cloud amounts
     &  , frac_cloud(nd_profile, nd_layer, nd_cloud_type)               &
!           Cloud amounts
     &  , condensed_min_dim(nd_cloud_component)                         &
!           Minimum dimensions of condensed components
     &  , condensed_max_dim(nd_cloud_component)
!           Maximum dimensions of condensed components
!
!     Properties of aerosols:
      logical                                                           &
     &    l_aerosol                                                     &
!           Flag to enable direct aerosol effects within
!           the radiation code
     &  , l_aerosol_ccn
!           Flag to define CCN from aerosols
      real                                                              &
     &    aerosol_mix_ratio(nd_profile, nd_layer                        &
     &      , lw_spectrum%npd_aerosol_species)
!           Mixing ratios of aerosols
!
!     Local algorithmic variables:
      real                                                              &
     &    euler_factor
!           Weighting applied to the last term of alternating series
      integer                                                           &
     &    i_scatter_method_band(lw_spectrum%npd_band)
!           Method of treating scattering in the calculation
!           of optical propeties in each band
!
!     Gathered viewing directions:
      real                                                              &
     &    viewing_direction_g(nd_radiance_profile, nd_direction, 2)
!           Satellite viewing directions
!
!     Fluxes:
      real                                                              &
     &    flux_direct(nd_flux_profile, 0: nd_layer, nd_channel)         &
!           Direct flux
     &  , flux_direct_clear(nd_flux_profile, 0: nd_layer, nd_channel)   &
!           Clear-sky direct flux
     &  , flux_net(nd_flux_profile, 0: nd_layer, nd_channel)            &
!           Downward/net flux
     &  , flux_net_clear(nd_flux_profile, 0: nd_layer, nd_channel)      &
!           Clear-sky downward/net flux
     &  , flux_up(nd_flux_profile, 0: nd_layer, nd_channel)             &
!           Upward flux
     &  , flux_up_clear(nd_flux_profile, 0: nd_layer, nd_channel)       &
!           Clear-sky upward flux
     &  , surface_down_flux_g(nd_profile)                               &
!           Downward flux at surface
     &  , surf_down_clr_g(nd_profile)
!           Clear-sky downward flux at surface
!
!     local arrays to fill diagnostics:
      real :: total_cloud_cover_g(nd_profile)
!             cloud fraction at gathered points
!     Tiled fluxes:
      real                                                              &
     &    flux_up_tile(nd_point_tile, nd_tile, nd_channel)
!           Upward fluxes at tiled surface points
!
!     Radiances:
      real                                                              &
     &    radiance(nd_radiance_profile, nd_viewing_level                &
     &      , nd_direction, nd_channel)
!           Radiances calculated
!
!     Fields required for call to radiation code but not used
      INTEGER, Parameter :: nd_j_profile=1
      INTEGER :: i_gas
!
!
!     Auxiliary variables:
      real                                                              &
     &    dacon                                                         &
!           Difference in A's
     &  , dbcon                                                         &
!           Difference in B's
     &  , weight_band(lw_spectrum%npd_band)                             &
!           Weighting factors for bands
     &  , nullmmr
!           Null mass mixing ratio
      parameter(nullmmr=0.0)
!
!     Dummy fields for radiation code
      logical                                                           &
     &     l_dummy
      real                                                              &
     &     dummy                                                        &
     &    ,dummy1d(1)                                                   &
     &    ,dummy2d(1,1)                                                 &
     &    ,dummy3d(1,1,1)
!
!
!     Subroutines called:
      external                                                          &
     &    r2_set_gas_mix_ratio, r2_set_thermodynamic                    &
     &  , r2_set_aerosol_field, r2_set_cloud_field                      &
     &  , r2_set_cloud_parametrization                                  &
     &  , r2_set_surface_field_lw, radiance_calc                        &
     &  , r2_calc_total_cloud_cover
!
!
!
!
!     Initialize the error flag for the radiation code.
      ierr=i_normal
#if defined(RAD_DBG)
! DEPENDS ON: r2_lw_dump
      call r2_lw_dump(ierr                                              &
!                       Physical dimensions
     &  , LW_spectrum%npd_band                                          &
     &  , n_points, nlevs, n_layer, nclds                               &
     &  , nwet, nozone, row_length, rows, nd_field                      &
     &  , nd_field_flux_diag, nd_field_rad_diag                         &
     &  , nd_profile, nd_layer, nd_column, n_cca_lev, nd_channel        &
     &  , nd_flux_profile, nd_radiance_profile                          &
     &  , nd_viewing_level, nd_direction                                &
     &  , nd_cloud_component, nd_cloud_type                             &
     &  , nd_brdf_basis_fnc, nd_brdf_trunc                              &
     &  , nd_point_tile, nd_tile, id_ct                                 &
!                       Gaseous mixing ratios
     &  , h2o, co2, o3                                                  &
     &  , co2_dim1, co2_dim2, co2_3d, l_co2_3d                          &
     &  , n2o_mix_ratio, ch4_mix_ratio                                  &
     &  , cfc11_mix_ratio, cfc12_mix_ratio, cfc113_mix_ratio            &
     &  , hcfc22_mix_ratio, hfc125_mix_ratio, hfc134a_mix_ratio         &
!                       Thermodynamic variables
     &  , tac,  tstar,  tstar_solid, tstar_sea, l_ctile, pstar          &
     &  , p_layer_boundaries                                            &
     &  , p_layer_centres                                               &
     &  , height_theta                                                  &
     &  , height_rho                                                    &
!                       Options for treating clouds
     &  , global_cloud_top                                              &
!                       Stratiform cloud fields
     &  , l_cloud_water_partition                                       &
     &  , lca_area, lca_bulk, lccwc1, lccwc2                            &
!                       Convective cloud fields
     &  , cca, cccwp, ccb, cct                                          &
!                       Surface fields
     &  , land, flandg, ice_fraction                                    &
     &  , lying_snow                                                    &
!                       Aerosol fields
     &  , l_climat_aerosol, l_clim_aero_hgt, bl_depth, n_levels_bl      &
     &  , l_use_sulpc_direct, l_use_sulpc_indirect                      &
     &  , sulp_dim1,sulp_dim2                                           &
     &  , accum_sulphate, aitken_sulphate, diss_sulphate                &
     &  , sea_salt_film, sea_salt_jet, l_use_seasalt_indirect           &
     &  , l_use_seasalt_direct, salt_dim_a, salt_dim_b                  &
     &  , l_use_soot_direct, soot_dim1, soot_dim2                       &
     &  , fresh_soot, aged_soot                                         &
     &  , aero_meso, l_murk_rad                                         &
!                       Level of tropopause
     &  , trindx                                                        &
!                       Algorithmic options
     &  , list                                                          &
     &  , pts, l_scale_inc                                              &
     &  , i_segment, i_call                                             &
!                       Satellite viewing geometry
     &  , n_viewing_direction, viewing_direction, n_viewing_level       &
     &  , viewing_level                                                 &
!                       Diagnostics
     &  , n_channel, map_channel                                        &
     &  , row_list, col_list                                            &
     &  )
#endif
!
!     Becuase the radiation code may be invoked on subsamples of
!     the same segment with different selections of points there
!     may be ocassions when this routine is called only to fill
!     diagnostics with zeros.
      if (n_points >  0) then
!
!       Set the logical flag for dummy diagnostics not available from
!       the lower code in the long-wave to .FALSE..
        l_dummy=.false.
!
!
!
!       Set dynamic array sizes and dependent options.
! DEPENDS ON: r2_set_option_lw
        call r2_set_option_lw(ierr                                      &
     &    , n_layer, lw_spectrum%n_band                                 &
     &    , lw_control%isolir                                           &
     &    , lw_control%i_gas_overlap, i_gas_overlap                     &
     &    , lw_control%l_aerosol, l_climat_aerosol                      &
     &    , l_use_sulpc_direct, l_use_soot_direct, l_use_biogenic       &
     &    , l_use_seasalt_direct, l_murk_rad, l_aerosol                 &
     &    , l_use_sulpc_indirect, l_aerosol_ccn, n_arcl_species         &
     &    , lw_control%l_global_cloud_top                               &
     &    , global_cloud_top, n_cloud_top_global                        &
     &    , lw_diag%l_clear_olr, lw_diag%l_surf_down_clr                &
     &    , lw_diag%l_clear_hr, l_clear                                 &
     &    , lw_control%i_angular_integration                            &
     &    , lw_control%i_solver, i_solver_clear                         &
     &    , lw_control%l_rescale, n_order_forward                       &
     &    , lw_control%i_truncation, lw_control%ls_global_trunc         &
     &    , lw_control%l_euler_trnf, euler_factor                       &
     &    , lw_control%i_scatter_method, i_scatter_method_band          &
     &    , l_rad_tile                                                  &
     &    , weight_band                                                 &
     &    , nd_overlap_coeff, nd_2sg_profile, nd_layer_clr              &
     &    , nd_source_coeff, nd_max_order, nd_sph_coeff                 &
     &    , nd_region                                                   &
     &    , nd_profile, lw_spectrum%npd_band                            &
     &    )
!
!       If radiances are to be calculated the viewing directions
!       are to be copied into the gathered array.
        if ( (lw_control%i_angular_integration ==                       &
     &        ip_spherical_harmonic).and.                               &
     &       (lw_control%i_sph_mode == ip_sph_mode_rad) ) then
          do ll=1, n_points
            l=list(ll)
            viewing_direction_g(ll, 1, 1)=viewing_direction1(l, 1)
            viewing_direction_g(ll, 1, 2)=viewing_direction2(l, 1)
          enddo
        endif
!
!
!
!       Set the mixing ratios of gases.
! DEPENDS ON: r2_set_gas_mix_ratio
        call r2_set_gas_mix_ratio(ierr, i_call                          &
     &    , n_points, nlevs, n_layer, nwet, nozone                      &
     &    , list, lw_control%l_extra_top                                &
     &    , lw_spectrum%n_absorb, lw_spectrum%type_absorb               &
     &    , lw_control%l_n2o, lw_control%l_ch4, lw_control%l_cfc11      &
     &    , lw_control%l_cfc12, .false.                                 &
     &    , lw_control%l_cfc113, lw_control%l_hcfc22                    &
     &    , lw_control%l_hfc125, lw_control%l_hfc134a                   &
     &    , h2o, co2, o3, n2o_mix_ratio, ch4_mix_ratio                  &
     &    , cfc11_mix_ratio, cfc12_mix_ratio, nullmmr                   &
     &    , cfc113_mix_ratio, hcfc22_mix_ratio, hfc125_mix_ratio        &
     &    , hfc134a_mix_ratio                                           &
     &    , gas_mix_ratio                                               &
     &    , co2_dim1, co2_dim2, co2_3d, l_co2_3d                        &
     &    , nd_field, nd_profile, nd_layer, lw_spectrum%npd_species     &
     &    )
        if (ierr /= i_normal) return
!
!
!       Calculate pressures and temperatures.
! DEPENDS ON: r2_set_thermodynamic
        call r2_set_thermodynamic(n_points,nlevs,n_layer,nwet,list      &
     &    , lw_control%l_extra_top, .true.                              &
     &    , pstar, tstar, tstar_solid, tstar_sea                        &
     &    , p_layer_boundaries                                          &
     &    , p_layer_centres                                             &
     &    , height_theta                                                &
     &    , height_rho                                                  &
     &    , tac                                                         &
     &    , rho_r2, r_rho_levels, r_theta_levels                        &
     &    , q, qcl, qcf, qcf2, qrain, qgraup                            &
     &    , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio       &
     &    , p, t, t_bdy, t_surface, t_solid, t_sea, d_mass              &
     &    , layer_heat_capacity                                         &
     &    , nd_field, nd_profile, nd_layer                              &
     &    )
!
!
!       set sea-salt array dimensions.
        if (l_use_seasalt_direct) then
           salt_dim_dir_a=salt_dim_a
           salt_dim_dir_b=salt_dim_b
        else
           salt_dim_dir_a=1
           salt_dim_dir_b=1
        endif
!
!
!       Set the mixing ratios of aerosols.
        if (l_aerosol) then
! DEPENDS ON: r2_set_aerosol_field
          call r2_set_aerosol_field(ierr                                &
     &      , n_points, nlevs, n_layer, lw_spectrum%n_aerosol           &
     &      , lw_spectrum%type_aerosol                                  &
     &      , list, lw_control%l_extra_top                              &
     &      , l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero    &
     &      , bl_depth, t, n_levels_bl, l_murk_rad, aero_meso           &
     &      , l_use_sulpc_direct                                        &
     &      , sulp_dim1, sulp_dim2                                      &
     &      , accum_sulphate, aitken_sulphate                           &
     &      , l_use_seasalt_direct, salt_dim_dir_a, salt_dim_dir_b      &
     &      , sea_salt_film, sea_salt_jet, p                            &
     &      , l_use_soot_direct, soot_dim1, soot_dim2                   &
     &      , fresh_soot, aged_soot                                     &
     &      , l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic    &
     &      , n_arcl_species, n_arcl_compnts, i_arcl_compnts            &
     &      , l_use_arcl, arcl_dim1, arcl_dim2, arcl                    &
     &      , land, lying_snow, pstar                                   &
     &      , p_layer_boundaries                                        &
     &      , trindx                                                    &
     &      , aerosol_mix_ratio                                         &
     &      , nd_field, nd_profile, nd_layer                            &
     &      , lw_spectrum%npd_aerosol_species                           &
     &      )
        endif
!
!
!       Assign the properties of clouds. A dummy array must be passed
!       for the microphysical diagnostics since they are not available
!       through STASH in the long-wave.
!
! DEPENDS ON: r2_set_cloud_parametrization
        call r2_set_cloud_parametrization(ierr, lw_spectrum%n_band      &
     &    , lw_control%i_st_water, lw_control%i_cnv_water               &
     &    , lw_control%i_st_ice, lw_control%i_cnv_ice                   &
     &    , lw_spectrum%l_drop_type                                     &
     &    , lw_spectrum%i_drop_parametrization                          &
     &    , lw_spectrum%n_drop_phf_term, lw_spectrum%drop_parameter_list&
     &    , lw_spectrum%drop_parm_min_dim, lw_spectrum%drop_parm_max_dim&
     &    , lw_spectrum%l_ice_type                                      &
     &    , lw_spectrum%i_ice_parametrization                           &
     &    , lw_spectrum%n_ice_phf_term,lw_spectrum%ice_parameter_list   &
     &    , lw_spectrum%ice_parm_min_dim, lw_spectrum%ice_parm_max_dim  &
     &    , i_condensed_param, condensed_n_phf, condensed_param_list    &
     &    , condensed_min_dim, condensed_max_dim                        &
     &    , lw_spectrum%npd_band, lw_spectrum%npd_drop_type             &
     &    , lw_spectrum%npd_ice_type, lw_spectrum%npd_cloud_parameter   &
     &    , nd_cloud_component                                          &
     &    )
        if (ierr /= i_normal) return
!
!
!       set sea-salt array dimensions.
        if (l_use_seasalt_indirect) then
           salt_dim_ind_a=salt_dim_a
           salt_dim_ind_b=salt_dim_b
        else
           salt_dim_ind_a=1
           salt_dim_ind_b=1
        endif
!
!
!       the gathered land flag is required for cloud microphysics.
!       It may be more logical to put this in r2_set_surface_field_lw
!       and move that earlier in the calling sequence.
        do l=1, n_points
          land_g(l)=land(list(l))
          flandg_g(l)=flandg(list(l))
        enddo
!
!
! DEPENDS ON: r2_set_cloud_field
        call r2_set_cloud_field(n_points, nlevs, n_layer, nclds         &
     &    , list                                                        &
     &    , p, t, d_mass                                                &
     &    , ccb, cct, cca, cccwp                                        &
     &    , lccwc1, lccwc2, lca_area, lca_bulk                          &
     &    , lw_control%l_microphysics, l_aerosol_ccn                    &
     &    , sea_salt_film, sea_salt_jet                                 &
     &    , l_use_seasalt_indirect, salt_dim_ind_a, salt_dim_ind_b      &
     &    , sulp_dim1, sulp_dim2, accum_sulphate, diss_sulphate         &
     &    , aitken_sulphate, lying_snow                                 &
     &    , l_cloud_water_partition,  land_g, flandg_g                  &
     &    , lw_control%i_cloud_representation, i_condensed_param        &
     &    , condensed_min_dim, condensed_max_dim                        &
     &    , n_condensed, type_condensed                                 &
     &    , w_cloud, n_cloud_type, frac_cloud                           &
     &    , lw_control%l_local_cnv_partition                            &
     &    , condensed_mix_ratio, condensed_dim_char                     &
!                       Microphysical diagnostics are not available
!                       in this spectral region.
     &    , dummy2d, .false., dummy2d, .false.                          &
     &    , dummy2d, .false., dummy2d, .false.                          &
     &    , dummy2d, .false.                                            &
     &    , dummy2d, .false., dummy2d, .false.                          &
     &    , dummy2d, .false., dummy2d, .false.                          &
     &    , col_list, row_list, row_length, rows                        &
     &    , nd_field, nd_profile, nd_layer, id_ct                       &
     &    , lw_spectrum%npd_aerosol_species                             &
     &    , nd_cloud_component, nd_cloud_type                           &
     &    , n_cca_lev, Ntot_land, Ntot_sea                              &
     &    )
!
!
! DEPENDS ON: r2_set_surface_field_lw
        call r2_set_surface_field_lw(ierr                               &
     &    , n_points, list                                              &
     &    , lw_spectrum%n_band, lw_control%ls_brdf_trunc                &
     &    , l_ctile, flandg                                             &
     &    , n_brdf_basis_fnc, f_brdf, rho_alb                           &
     &    , land, ice_fraction, tstar_sea, tstar_solid                  &
     &    , l_rad_tile, n_point_tile, n_tile, list_tile, rho_alb_tile   &
     &    , frac_tile, t_tile, list_tile_outer, index_tile              &
     &    , nd_field, nd_profile, lw_spectrum%npd_band                  &
     &    , nd_brdf_basis_fnc                                           &
     &    , nd_brdf_trunc, nd_point_tile, nd_tile                       &
     &    )
!
!
!
        radiance = 0.0
!
! DEPENDS ON: radiance_calc
        call radiance_calc(ierr                                         &
!                       Logical flags for processes
     &    , lw_control%l_rayleigh, l_aerosol                            &
     &    , lw_control%l_gas, lw_control%l_continuum                    &
     &    , lw_control%l_cloud, lw_control%l_drop, lw_control%l_ice     &
!                       Angular integration
     &    , lw_control%i_angular_integration, lw_control%l_rescale      &
     &    , n_order_forward                                             &
     &    , lw_control%i_2stream, lw_control%n_order_gauss              &
     &    , lw_control%i_truncation, lw_control%ls_global_trunc         &
     &    , lw_control%ms_min, lw_control%ms_max                        &
     &    , lw_control%accuracy_adaptive, euler_factor                  &
     &    , lw_control%l_henyey_greenstein_pf, lw_control%ls_brdf_trunc &
     &    , lw_control%i_sph_algorithm, lw_control%n_order_phase_solar  &
     &    , n_viewing_direction, viewing_direction_g                    &
     &    , n_viewing_level, viewing_level                              &
     &    , lw_control%i_sph_mode                                       &
!                       Treatment of scattering
     &    , i_scatter_method_band                                       &
!                       Options for treating clouds
     &    , lw_control%l_global_cloud_top, n_cloud_top_global           &
     &    , l_inhom_cloud, inhom_cloud                                  &
!                       Options for solver
     &    , lw_control%i_solver                                         &
!                       Properties of diagnostics
     &    , map_channel                                                 &
!                       General spectral properties
     &    , lw_spectrum%n_band                                          &
     &    , lw_control%first_band, lw_control%last_band                 &
     &    , weight_band                                                 &
!                       General atmospheric properties
     &    , n_points, n_layer                                           &
     &    , p, t, t_surface, t_bdy, d_mass                              &
!                       Spectral region
     &    , lw_control%isolir                                           &
!                       Solar fields
     &    , dummy, dummy, lw_spectrum%solar_flux_band                   &
     &    , lw_spectrum%rayleigh_coefficient                            &
!                       Infra-red fields
     &    , lw_spectrum%n_deg_fit, lw_spectrum%thermal_coefficient      &
     &    , lw_spectrum%t_ref_planck, lw_control%l_ir_source_quad       &
!                       Gaseous absorption
     &    , i_gas_overlap, i_gas, gas_mix_ratio                         &
     &    , lw_spectrum%n_band_absorb, lw_spectrum%index_absorb         &
     &    , lw_spectrum%i_band_esft                                     &
     &    , lw_spectrum%w_esft, lw_spectrum%k_esft                      &
     &    , lw_spectrum%i_scale_esft, lw_spectrum%i_scale_fnc           &
     &    , lw_spectrum%scale_vector                                    &
     &    , lw_spectrum%p_reference, lw_spectrum%t_reference            &
     &    , l_mod_k_flux                                                &
!                       Doppler broadening
     &    , lw_spectrum%l_doppler_present                               &
     &    , lw_spectrum%doppler_correction                              &
!                       Surface fields
     &    , n_brdf_basis_fnc, rho_alb, f_brdf                           &
!                       Tiling options for heterogeneous surfaces
     &    , l_rad_tile, n_point_tile, n_tile, list_tile, rho_alb_tile   &
     &    , frac_tile, t_tile                                           &
!                       Continuum absorption
     &    , lw_spectrum%n_band_continuum                                &
     &    , lw_spectrum%index_continuum, lw_spectrum%index_water        &
     &    , lw_spectrum%k_continuum, lw_spectrum%i_scale_fnc_cont       &
     &    , lw_spectrum%scale_continuum                                 &
     &    , lw_spectrum%p_ref_continuum, lw_spectrum%t_ref_continuum    &
!                       Properties of aerosols
     &    , lw_spectrum%n_aerosol, aerosol_mix_ratio                    &
     &    , lw_spectrum%aerosol_absorption                              &
     &    , lw_spectrum%aerosol_scattering                              &
     &    , lw_spectrum%n_aerosol_phf_term                              &
     &    , lw_spectrum%aerosol_phase_fnc                               &
     &    , lw_spectrum%i_aerosol_parametrization                       &
     &    , lw_spectrum%nhumidity, lw_spectrum%humidities               &
!                       Properties of clouds
     &    , n_condensed, type_condensed                                 &
     &    , lw_control%i_cloud, lw_control%i_cloud_representation       &
     &    , w_cloud, n_cloud_type, frac_cloud                           &
     &    , condensed_mix_ratio, condensed_dim_char                     &
     &    , i_condensed_param, condensed_n_phf, condensed_param_list    &
     &    , dp_corr_strat, dp_corr_conv                                 &
!                       Calculated Fluxes
     &    , flux_direct, dummy3d, flux_net, flux_up                     &
     &    , dummy3d, dummy3d, dummy3d                                   &
     &    , l_dummy, l_dummy, l_dummy, l_dummy                          &
!                       Calculated Radiances
     &    , radiance, dummy2D                                           &
!                       Options for clear-sky fluxes
     &    , l_clear, i_solver_clear                                     &
!                       Clear-sky fluxes calculated
     &    , flux_direct_clear, flux_net_clear, flux_up_clear            &
!                       Special Surface Fluxes
     &    , .false., dummy1d, dummy1d, dummy1d, dummy1d, dummy1d        &
!                       Tiled Surface Fluxes
     &    , flux_up_tile, dummy3d                                       &
!                       Dimensions of arrays
     &    , nd_profile, nd_layer, nd_column, nd_layer_clr, id_ct        &
     &    , nd_2sg_profile, nd_flux_profile, nd_radiance_profile        &
     &    , nd_j_profile, nd_channel                                    &
     &    , lw_spectrum%npd_band, lw_spectrum%npd_species               &
     &    , lw_spectrum%npd_esft_term, lw_spectrum%npd_scale_variable   &
     &    , lw_spectrum%npd_continuum                                   &
     &    , lw_spectrum%npd_aerosol_species, lw_spectrum%npd_humidities &
     &    , lw_spectrum%npd_cloud_parameter                             &
     &    , lw_spectrum%npd_thermal_coeff, nd_source_coeff              &
     &    , nd_brdf_basis_fnc, nd_brdf_trunc                            &
     &    , lw_spectrum%npd_phase_term, nd_max_order, nd_sph_coeff      &
     &    , nd_direction, nd_viewing_level                              &
     &    , nd_region, nd_cloud_type, nd_cloud_component                &
     &    , nd_overlap_coeff, nd_point_tile, nd_tile                    &
!                       Number of call
     &    , i_call                                                      &
     &    )
        if (ierr /= i_normal) return
!
!
      endif
!
!     Assignment of diagnostics:
!
!     Processing depends on whether the code has been invoked to
!     calculate radiances or fluxes.
      if ( (lw_control%i_angular_integration == ip_two_stream).or.      &
     &     ( (lw_control%i_angular_integration ==                       &
     &        ip_spherical_harmonic).and.                               &
     &        (lw_control%i_sph_mode == ip_sph_mode_flux) ) ) then
!
!       Convert downward fluxes to net fluxes.
        do i=0, n_layer
          do l=1, n_points
            flux_net(l, i, 1)=flux_net(l, i, 1)-flux_up(l, i, 1)
          enddo
        enddo
        if (l_clear) then
          do i=0, n_layer
            do l=1, n_points
              flux_net_clear(l, i, 1)                                   &
     &          =flux_net_clear(l, i, 1)-flux_up_clear(l, i, 1)
            enddo
          enddo
        endif
!
!
!       OLR:
!
        do l=1, n_points
          olr(list(l))=-flux_net(l, 0, 1)
        enddo
        if (LW_diag%L_clear_olr) then
          do l=1, n_points
            LW_diag%clear_olr(col_list(l), row_list(l))                 &
     &        =-flux_net_clear(l, 0, 1)
          enddo
        endif
!
!
!       Total cloud cover:
!
        if (LW_diag%l_total_cloud_cover) then
          write(iu_err, *)                                              &
     &      '*** warning total cloud cover is not yet fully defined.'
!! DEPENDS ON: r2_calc_total_cloud_cover
!          call r2_calc_total_cloud_cover(n_points, n_layer, nclds       &
!     &      , lw_control%i_cloud, w_cloud, total_cloud_cover_g          &
!     &      , nd_profile, nd_layer, id_ct                               &
!     &      )
          Do l=1, n_points
            LW_diag%total_cloud_cover(col_list(l), row_list(l))         &
!     &        =total_cloud_cover_g(l)
     &        =0.0
          Enddo
        endif
!
!
!       Net flux at the tropopause:
!
        if (LW_diag%l_net_flux_trop) then
          do l=1, n_points
            LW_diag%net_flux_trop(col_list(l), row_list(l))             &
     &        =flux_net(l, n_layer+1-trindx(l), 1)
          enddo
        endif
!
!
!       Downward flux at the tropopause:
!
        if (LW_diag%l_down_flux_trop) then
          do l=1, n_points
            LW_diag%down_flux_trop(col_list(l), row_list(l))            &
     &        =flux_net(l, n_layer+1-trindx(l), 1)                      &
     &        +flux_up(l, n_layer+1-trindx(l), 1)
          enddo
        endif
!
!
!       Downward flux at the surface:
!
        if (LW_diag%l_surface_down_flux) then
          do l=1, n_points
            LW_diag%surface_down_flux(col_list(l), row_list(l))         &
     &        =flux_net(l, n_layer, 1)                                  &
     &        +flux_up(l, n_layer, 1)
          enddo
        endif
!
!
!       Clear-sky downward flux at the surface:
!
        if (LW_diag%l_surf_down_clr) then
          do l=1, n_points
            LW_diag%surf_down_clr(col_list(l), row_list(l))             &
     &        =flux_net_clear(l, n_layer, 1)                            &
     &        +flux_up_clear(l, n_layer, 1)
          enddo
        endif
        
!
! ######################################################
! CLOUD WATER MIXING RATIOS
! ######################################################
!
!     LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (LIQUID)
!     ==============================================================
      If (LW_diag%L_ls_qcl_rad) Then
        Do I=1, N_LAYER+1-id_ct
          Do L=1, n_points
            LW_diag%ls_qcl_rad(col_list(L),row_list(L),I)               &
     &         = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_ST_WATER)   &
     &         * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SW)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (ICE)
!     ==============================================================
      If (LW_diag%L_ls_qcf_rad) Then
        Do I=1, N_LAYER+1-id_ct
          Do L=1, n_points
            LW_diag%ls_qcf_rad(col_list(L),row_list(L),I)               &
     &         = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_ST_ICE)     &
     &         * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SI)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (LIQUID)
!     ==============================================================
      If (LW_diag%L_cc_qcl_rad) Then
        Do I=1, N_LAYER+1-id_ct
          Do L=1, n_points
            LW_diag%cc_qcl_rad(col_list(L),row_list(L),I)               &
     &         = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_CNV_WATER)  &
     &         * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CW)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (ICE)
!     ==============================================================

      If (LW_diag%L_cc_qcf_rad) Then
        Do I=1, N_LAYER+1-id_ct
          Do L=1, n_points
            LW_diag%cc_qcf_rad(col_list(L),row_list(L),I)               &
     &         = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_CNV_ICE)    &
     &         * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CI)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
! ######################################################
! CLOUD FRACTIONS
! ######################################################
!
!     LARGE-SCALE cloud GRIDBOX FRACTION seen by radiation. (LIQUID)
!     ==============================================================
      If (LW_diag%L_ls_cl_rad) Then
        Do I=1, N_LAYER+1-id_ct
          Do L=1, n_points
            LW_diag%ls_cl_rad(col_list(L),row_list(L),I)                &
     &         = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SW)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     LARGE-SCALE cloud GRIDBOX fraction seen by radiation. (ICE)
!     ==============================================================
      If (LW_diag%L_ls_cf_rad) Then
        Do I=1, N_LAYER+1-id_ct
          Do L=1, n_points
            LW_diag%ls_cf_rad(col_list(L),row_list(L),I)                &
     &         = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SI)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     CONVECTIVE cloud GRIDBOX fraction seen by radiation. (LIQUID)
!     ==============================================================
      If (LW_diag%L_cc_cl_rad) Then
        Do I=1, N_LAYER+1-id_ct
          Do L=1, n_points
            LW_diag%cc_cl_rad(col_list(L),row_list(L),I)                &
     &         = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CW)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif
!
!
!     CONVECTIVE cloud GRIDBOX FRACTION seen by radiation. (ICE)
!     ==============================================================
      If (LW_diag%L_cc_cf_rad) Then
        Do I=1, N_LAYER+1-id_ct
          Do L=1, n_points
            LW_diag%cc_cf_rad(col_list(L),row_list(L), I)               &
     &         = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CI)             &
     &         * W_CLOUD(L,N_LAYER+1-I)
          Enddo
        Enddo
      Endif     
!
!
!
!
!       Output arrays:
!
!       Convert the fluxes to increments in the heating rate except at
!       the surface: there, the net downward flux is assigned to LWOUT.
        do k=nlevs, 1, -1
          if (l_scale_inc) then
            do l=1, n_points
              lwout(list(l), k+1)=(flux_net(l, n_layer-k, 1)            &
     &          -flux_net(l, n_layer+1-k, 1))                           &
     &          *pts/layer_heat_capacity(l, n_layer+1-k)
            enddo
          else
            do l=1, n_points
              lwout(list(l), k+1)=(flux_net(l, n_layer-k, 1)            &
     &          -flux_net(l, n_layer+1-k, 1))                           &
     &          /layer_heat_capacity(l, n_layer+1-k)
            enddo
          endif
          if (LW_diag%l_clear_hr) then
!           The factor of PTS is included here to yield a rate from an
!           increment.
            do l=1, n_points
              LW_diag%clear_hr(col_list(l), row_list(l), k)             &
     &          =(flux_net_clear(l, n_layer-k, 1)                       &
     &          -flux_net_clear(l, n_layer+1-k, 1))                     &
     &          /layer_heat_capacity(l, n_layer+1-k)
            enddo
          endif
        enddo
!
        if (lw_control%l_extra_top) then
!         calculate the radiation absorbed in the extra layer
!         above the top of the rest of the model.
          do l=1, n_points
            top_absorption(list(l))=flux_net(l, 0, 1)                   &
     &        -flux_net(l, n_layer-nlevs, 1)
          enddo
        endif
!
        do l=1, n_points
          lwout(list(l), 1)=flux_net(l, n_layer, 1)
        enddo
!
!       Separate the contributions over open sea and sea-ice.
!       Fluxes returned from the radiation code itself are not
!       weighted by the fraction of the tile, but here are converted
!       to grid-box mean values. This split is possible only if the
!       ocean surface has been tiled.
!
        if (l_rad_tile) then
!
!         The variable flandg is set even if coastal tiling is not
!         used, so fairly generic code can be written.
!
!         It is simplest to zero LWsea at all points and reset
!         elsewhere.
          lwsea(list(1:n_points)) = 0.0
!
          do ll=1, n_points
            l=list(ll)
            if ( (flandg(l) < TINY(flandg)) .AND.                       &
     &           (ice_fraction(l) < TINY(ice_fraction) )                &
     &         ) then
!             This point is open sea with no sea ice.
              lwsea(l)=lwout(l, 1)
              lwout(l, 1)=0.0
            endif
          enddo
!
!         Tiled points will have both land and sea. Note that the
!         channel index of flux_up_tile is hard-wired to 1 because
!         we don't envisage calling the code in other cases.
          do lll=1, n_point_tile
            ll=list_tile(lll)
            l=list_tile_outer(lll)
            lwsea(l)=(1.0-ice_fraction(l))*(flux_net(ll, n_layer, 1)    &
     &        +flux_up(ll, n_layer, 1)                                  &
     &        -flux_up_tile(lll, index_tile(ip_ocean_tile), 1))
            lwout(l, 1)=lwout(l, 1)-(1.0-flandg(l))*lwsea(l)
          enddo
!
!         The remaining points are entirely land points and lwout
!         need not be altered.
!
        else
!
!         Without radiative tiling we must assume that fluxes are
!         uniform across the grid-box.
          Where (flandg(list(1:n_points)) < 1.0-TINY(flandg))
            lwsea(list(1:n_points))=(1.0-ice_fraction(list(1:n_points)))&
     &        *flux_net(1:n_points, n_layer, 1)
            lwout(list(1:n_points), 1)=lwout(list(1:n_points), 1)       &
     &        -(1.0-flandg(list(1:n_points)))*lwsea(list(1:n_points))
          Endwhere
!
        endif
!
      else
!
!       Radiances are being calculated.
!
        if (lw_diag%l_toa_radiance) then
          do ic=1, n_channel
            do l=1, n_points
              lw_diag%toa_radiance(col_list(l), row_list(l), ic)        &
     &          =radiance(l, 1, 1, ic)
            enddo
          enddo
        endif
!
      endif
!
!
!
      return
      END SUBROUTINE r2_lwrad3z
#endif
