#if defined(A01_3Z) && defined(RAD_DBG)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to dump inputs to shortwave radiation.
!
! Method:
!       The inputs to the shortwave radiation are written to a binary
!       file whence they may later be read for debugging.
!
! Current owner of code: J.-C. Thelen
!
! Description of code:
!   Fortran 90
!
!- ---------------------------------------------------------------------
      subroutine r2_sw_dump(ierr                                        &
!                     physical dimensions
     &  , nd_band                                                       &
     &  , nlit, n_points, nlevs, n_layer, nclds                         &
     &  , nwet, nozone, row_length, rows                                &
     &  , nd_field, nd_field_flux_diag, nd_field_rad_diag               &
     &  , nd_profile, nd_layer, nd_column                               &
     &  , n_cca_lev, nd_channel, nd_flux_profile                        &
     &  , nd_radiance_profile, nd_viewing_level, nd_direction           &
     &  , nd_cloud_component, nd_cloud_type                             &
     &  , nd_brdf_basis_fnc, nd_brdf_trunc                              &
     &  , nd_point_tile, nd_tile, id_ct                                 &
!                       mixing ratios
     &  , h2o, co2, o3, o2_mix_ratio                                    &
     &  , co2_dim1, co2_dim2, co2_3d, l_co2_3d                          &
!                       pressure fields
     &  , pstar                                                         &
     &  , p_layer_boundaries                                            &
     &  , p_layer_centres                                               &
!                       temperatures
     &  , tac                                                           &
!                       options for treating clouds
     &  , global_cloud_top, l_microphysics                              &
!                       stratiform cloud fields
     &  , l_cloud_water_partition                                       &
     &  , lca_area, lca_bulk, lccwc1, lccwc2                            &
!                       convective cloud fields
     &  , cca, cccwp, ccb, cct                                          &
!                       surface fields
     &  , land_albedo, l_moses_ii, l_ctile                              &
     &  , land_alb, sice_alb, flandg                                    &
     &  , open_sea_albedo, ice_fraction, land, land0p5, lying_snow      &
!                       solar fields
     &  , coszin, lit, list, scs                                        &
!                       aerosol fields
     &  , l_climat_aerosol, l_clim_aero_hgt, bl_depth, n_levels_bl      &
     &  , l_use_sulpc_direct, l_use_sulpc_indirect                      &
     &  , sulp_dim1, sulp_dim2                                          &
     &  , accum_sulphate, aitken_sulphate, diss_sulphate                &
     &  , sea_salt_film, sea_salt_jet, l_use_seasalt_indirect           &
     &  , l_use_seasalt_direct, salt_dim_a, salt_dim_b                  &
     &  , l_use_soot_direct, soot_dim1, soot_dim2                       &
     &  , fresh_soot, aged_soot                                         &
     &  , aero_meso, l_murk_rad                                         &
!                       Level of tropopause
     &  , trindx                                                        &
!                       Algorithmic options
     &  , pts, l_scale_inc                                              &
     &  , i_segment, i_call                                             &
!                       Satellite viewing geometry
     &  , n_viewing_direction, viewing_direction, n_viewing_level       &
     &  , viewing_level, n_channel, map_channel                         &
!                       diagnostics
     &  , row_list, col_list, l_flux_below_690nm_surf                   &
     &  )
!
!
!
      Use cv_cntl_mod, Only:                                            &
          lcv_3d_cca

      implicit none
!
!
!
!     dummy arguments
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
     &    nd_band                                                       &
!           Number of spectral bands
     &  , nd_field                                                      &
!           Field size in calling program
     &  , nd_field_flux_diag                                            &
!           Field size for flux diagnostics
     &  , nd_field_rad_diag                                             &
!           Field size for radiance diagnostics
     &  , nd_profile                                                    &
!           Size of array of profiles
     &  , nd_layer                                                      &
!           Array sizes for layers
     &  , nd_column
!           Number of columns per point
      integer                                                           &
                !, intent(in)
     &    nd_channel                                                    &
!           Size allocated for diagnostic spectral bands
     &  , nd_flux_profile                                               &
!           Size allocated for output flux profiles
     &  , nd_radiance_profile                                           &
!           Number of profiles where radiances are required
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
     &  , nd_tile                                                       &
!           Size allocated for surface tiles
     &  , id_ct
!           Top level in arrays of cloud properties
!
!     Actual sizes used:
      integer                                                           &
                !, intent(in)
     &    n_points                                                      &
!           Number of points to be diagnosed including unlit points
     &  , nwet                                                          &
!           Number of wet levels
     &  , nozone                                                        &
!           Number of levels with ozone
     &  , nlevs                                                         &
!           number of layers in the main model
     &  , n_layer                                                       &
!           number of layers seen in the radiation scheme
     &  , nclds                                                         &
!           Number of cloudy levels
     &  , n_levels_bl                                                   &
!           Number of layers occupied by boundary-layer aerosol
!           if l_clim_aero_hgt is false.
     &  , n_cca_lev
!           Number of convective cloud levels
!
      logical                                                           &
     &    l_scale_inc
!           Flag for scaling of heating rates to increments
!
      INTEGER :: i_segment
!       Segment in call
      INTEGER :: i_call
!       Number of call
!
!
!
!     Gaseous mixing ratios
      real                                                              &
                !, intent(in)
     &    h2o(nd_field, nwet)                                           &
!           Mass mixing ratio of water
     &  , co2                                                           &
!           Mass mixing ratio of co2
     &  , o3(nd_field, nozone)                                          &
!           Mass mixing ratios of ozone
     &  , o2_mix_ratio
!           Mass mixing ratio of oxygen
!
!     General atmospheric properties:
      real                                                              &
                !, intent(in)
     &    pstar(nd_field)                                               &
!           Surface pressures
     &  , p_layer_boundaries(nd_field,0:nlevs)                          &
!            pressure at boundaries of layers
     &  , p_layer_centres(nd_field,0:nlevs)                             &
!            pressure at centres of layers
     &  , tac(nd_field, nlevs)
!           Temperatures at centres of layers
!
!     Incident solar radiation:
      integer                                                           &
                !, intent(in)
     &    nlit                                                          &
!           Number of lit points
     &  , list(nd_field)
!           List of lit points
      real                                                              &
                !, intent(in)
     &    coszin(nd_field)                                              &
!           Cosines of zenith angle
     &  , scs                                                           &
!           Scaling of solar incident field
     &  , lit(nd_field)
!           Fraction of time point is lit
!
!     Microphysical flag:
      logical                                                           &
                !, intent(in)
     &    l_microphysics
!           Flag for parametrized microphysics
!
!     Options for treating clouds
      integer                                                           &
                !, intent(in)
     &    global_cloud_top
!           Global topmost cloudy layer
!
!     Properties of stratiform clouds:
      logical                                                           &
                !, intent(in)
     &    l_cloud_water_partition
!           Flag to use prognostic cloud ice contents
      real                                                              &
                !, intent(in)
     &    lccwc1(nd_field, nclds+1/(nclds+1))                           &
!           Nominal liquid water contents
     &  , lccwc2(nd_field, nclds+1/(nclds+1))                           &
!           Nominal ice water contents
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
!           Fraction of convective cloud

!     Aerosols:
      logical                                                           &
                !, intent(in)
     &    l_climat_aerosol                                              &
!           Flag for climatological aerosol
     &  , l_clim_aero_hgt                                               &
!           flag to use the depth of the boundary layer to set
!           the climatological aerosol
     &  , l_murk_rad
!           flag for mesoscale model aerosol
      logical                                                           &
                !, intent(in)
     &    l_use_sulpc_direct                                            &
!           Flag to use sulphur cycle for direct effect
     &  , l_use_sulpc_indirect                                          &
!           Flag to use sulphur cycle for indirect effect
     &  , l_use_soot_direct                                             &
!           Flag to include the direct effect of soot
     &  , l_use_seasalt_indirect                                        &
!           flag to use sea-salt for indirect effect
     &  , l_use_seasalt_direct
!           flag to use sea-salt for direct effect
      integer                                                           &
                !, intent(in)
     &    sulp_dim1,sulp_dim2                                           &
!           Dimensions for _sulphate arrays, (P_FIELD,P_LEVELS or 1,1)
     &  , soot_dim1, soot_dim2                                          &
!           Dimensions for soot arrays (P_FIELD,P_LEVELS or 1,1)
     &  , salt_dim_a, salt_dim_b
!           dimensions for salt arrays on input (salt_dim_a=p_field
!           and salt_dim_b=p_levels, or else 1,1)
      real                                                              &
                !, intent(in)
     &    accum_sulphate(sulp_dim1, sulp_dim2)                          &
!           Mass mixing ratio of accumulation mode aerosol
     &  , aitken_sulphate(sulp_dim1, sulp_dim2)                         &
!           Mass mixing ratio of aitken mode aerosol
     &  , diss_sulphate(sulp_dim1, sulp_dim2)                           &
!           Mixing ratio of dissolved sulphate
     &  , sea_salt_film(salt_dim_a, salt_dim_b)                         &
!             number concentration of film-mode sea-salt aerosol
     &  , sea_salt_jet(salt_dim_a, salt_dim_b)                          &
!             number concentration of jet-mode sea-salt aerosol
     &  , fresh_soot(soot_dim1, soot_dim2)                              &
!           Soot mixing ratios
     &  , aged_soot(soot_dim1, soot_dim2)                               &
!           Soot mixing ratios
     &  , bl_depth(nd_field)                                            &
!           depth of the boundary layer
     &  , aero_meso(nd_field, nlevs)
!           mixing ratio of 'urban' aerosol of mesoscale model
!
!     Carbon cycle:
      logical                                                           &
     &    l_co2_3d
!           Controls use of 3D CO2 field
      integer                                                           &
                !, intent(in)
     &    co2_dim1, co2_dim2
!           Dimensions for CO2 array, (P_FIELD,P_LEVELS or 1,1)
      real                                                              &
                !, intent(in)
     &    co2_3d(co2_dim1, co2_dim2)
!           Mass mixing ratio of carbon dioxide
!     Properties of the surface:
      Logical, Intent(IN) :: Land(nd_field)
!                             Land mask
      Logical, Intent(IN) :: Land0p5(nd_field)
!                             Land mask (TRUE if land fraction > 0.5)
      Logical, Intent(IN) :: L_flux_below_690nm_surf
!                             Flag to calculate flux at wavelengths
!                             shorter than 690 nm at the surface:
!                             This may be required as a diagnostic
!                             or for MOSES
      Logical, Intent(IN) :: L_Moses_II
!                             Surface fluxes are required for MOSES II
      Logical, Intent(IN) :: L_Ctile
!                             Switch for coastal tiling
      real                                                              &
                !, intent(in)
     &    ice_fraction(nd_field)                                        &
!             fraction of sea ice in sea portion of grid box
     &  , land_albedo(nd_field,4)                                       &
!             MOSES ii land surface albedo fields
!             (*,1) - direct beam visible
!             (*,2) - diffuse visible
!             (*,3) - direct beam near-ir
!             (*,4) - diffuse near-ir
     &  , land_alb(nd_field)                                            &
!             surface albedo of land
     &  , sice_alb(nd_field)                                            &
!             surface albedo of sea-ice
     &  , flandg(nd_field)                                              &
!             land fraction in grid box
     &  , open_sea_albedo(nd_field, 2)                                  &
!           Surface albedo field of open sea
!             (direct and diffuse components)
     &  , lying_snow(nd_field)
!           Mass loading of lying snow
!
!                       level of tropopause
      integer                                                           &
     &    trindx(nd_field)
!           The layer boundary of the tropopause
!
!     Increment of time:
      real                                                              &
                !, intent(in)
     &    pts
!           Time increment
!
!     Satellite viewing geometry
      INTEGER, Intent(IN) :: n_channel
!           Number of channels calculated simultaneously in one
!           call to the radiation code
      INTEGER, Intent(IN) :: map_channel(nd_band)
!           Mapping of bands in the spectral file to channels in the
!           diagnostic output
      INTEGER, Intent(IN) :: n_viewing_direction
!           Number of viewing directions
      INTEGER, Intent(IN) :: n_viewing_level
!           Number of levels where the radiance is calculated
      REAL, Intent(IN) :: viewing_direction(nd_field, nd_direction, 2)
!           Satellite viewing directions
      REAL, Intent(IN) :: viewing_level(nd_viewing_level)
!           Levels where radiances are calculated
!
      integer, intent(in) :: row_list(nd_field)
!                              list of row indices of lit points
      integer, intent(in) :: col_list(nd_field)
!                              list of column indices of lit points
!
!
!     local variables.
      CHARACTER  (LEN=80) ::  file_dump_base
!             base of name of dump-file
      CHARACTER  (LEN=80) ::  file_dump
!             name of dump-file
      INTEGER :: i
!             loop variable
      INTEGER :: k
!             loop variable
      INTEGER :: l
!             loop variable
      INTEGER :: j
!             loop variable
      INTEGER :: n_entry = 0
!             number of entries into routine
      INTEGER :: ic_dump_target = 1
!             target for dumping
      INTEGER, Parameter :: iu_dump = 76
!             unit for dumping
      INTEGER :: irec
!             record
      INTEGER :: length_base
!             length of base name
      INTEGER :: node
!             number of processor in use
      INTEGER :: node_h
!             hundreds in number of processor in use
      INTEGER :: node_t
!             tens in number of processor in use
      INTEGER :: node_u
!             units in number of processor in use
      INTEGER :: i_segment_h
!             hundreds in ordinal of segment
      INTEGER :: i_segment_t
!             tens in ordinal of segment
      INTEGER :: i_segment_u
!             units in ordinal of segment
      INTEGER :: i_call_h
!             hundreds in ordinal of call
      INTEGER :: i_call_t
!             tens in ordinal of call
      INTEGER :: i_call_u
!             units in ordinal of call
      INTEGER :: i_zero
!             position of zero in the collating sequence
!
      INTEGER :: my_pe
!
!
      save n_entry
!
!
!
!
      file_dump_base='/u/m20/data/meg/t20db/rdump/r1_s'
!
!     increment the entry count and dump if the target is hit.
      n_entry=n_entry+1
      if (n_entry > ic_dump_target) return
!
!     check the number of the segment.
      if (i_segment >  999) then
!
         print *, '*** error: too many segments '                       &
     &      //'for naming convention.'
         return
!
      endif
!
!     find the processor in use.
      node=my_pe()
!     this code allows for only 999 processors.
      if (node >  999) then
!
         print *, '*** error: too many processors '                     &
     &      //'for naming convention.'
         return
!
      else
!
!        build a name for the output file.
         length_base=80
         do while ( (file_dump_base(length_base: length_base) == ' ')   &
     &      .and.(length_base >  2) )
            length_base=length_base-1
         enddo
!
         node_h=node/100
         node_t=(node-100*node_h)/10
         node_u=node-100*node_h-10*node_t
!
         i_segment_h=i_segment/100
         i_segment_t=(i_segment-100*i_segment_h)/10
         i_segment_u=i_segment-100*i_segment_h-10*i_segment_t
!
         i_call_h=i_call/100
         i_call_t=(i_call-100*i_call_h)/10
         i_call_u=i_call-100*i_call_h-10*i_call_t
!
         i_zero=ichar('0')
         file_dump(1: length_base+12)=file_dump_base(1: length_base)    &
     &      //'_'//char(node_h+i_zero)                                  &
     &      //char(node_t+i_zero)//char(node_u+i_zero)                  &
     &      //'_'//char(i_segment_h+i_zero)                             &
     &      //char(i_segment_t+i_zero)//char(i_segment_u+i_zero)        &
     &      //'_'//char(i_call_h+i_zero)                                &
     &      //char(i_call_t+i_zero)//char(i_call_u+i_zero)
!
      endif
!
!
      open(unit=iu_dump, file=file_dump(1: length_base+12)              &
     &      , form='unformatted', recl=512, access='direct')
!
!
      irec=1
!
!     dimensions of arrays:
      write(iu_dump, rec=irec)                                          &
     &    nlit, n_points, nlevs, n_layer, nclds                         &
     &  , nwet, nozone, row_length, rows                                &
     &  , nd_field, nd_field_flux_diag, nd_field_rad_diag               &
     &  , nd_profile, nd_layer, nd_column                               &
     &  , n_cca_lev, nd_channel, nd_flux_profile                        &
     &  , nd_radiance_profile, nd_viewing_level, nd_direction           &
     &  , nd_cloud_component, nd_cloud_type                             &
     &  , nd_brdf_basis_fnc, nd_brdf_trunc                              &
     &  , nd_point_tile, nd_tile, id_ct

      irec=irec+1
!
!     gaseous mixing ratios:
      do i=1, nwet
         do k=1, n_points/64
            write(iu_dump, rec=irec) (h2o(l, i)                         &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (h2o(l, i)                         &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
      write(iu_dump, rec=irec) co2
      irec=irec+1
      do i=1, nozone
         do k=1, n_points/64
            write(iu_dump, rec=irec) (o3(l, i), l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (o3(l, i)                          &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
      write(iu_dump, rec=irec) o2_mix_ratio
      irec=irec+1
      write(iu_dump, rec=irec) co2_dim1, co2_dim2, l_co2_3d
      irec=irec+1
      do i=1, co2_dim2
         do k=1, co2_dim1/64
            write(iu_dump, rec=irec) (co2_3d(l, i)                      &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (co2_dim1 >  64*(co2_dim1/64))  then
            write(iu_dump, rec=irec) (co2_3d(l, i)                      &
     &            , l=64*(co2_dim1/64)+1, co2_dim1)
            irec=irec+1
         endif
      enddo
!
!     thermodynamic fields:
      do k=1, n_points/64
         write(iu_dump, rec=irec) (pstar(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (pstar(l), l=64*(n_points/64)+1       &
     &         , n_points)
         irec=irec+1
      endif
      do i=0, nlevs
         do k=1, n_points/64
            write(iu_dump, rec=irec) (p_layer_boundaries(l, i)          &
     &        , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (p_layer_boundaries(l, i)          &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
      do i=0, nlevs
         do k=1, n_points/64
            write(iu_dump, rec=irec) (p_layer_centres(l, i)             &
     &        , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (p_layer_centres(l, i)             &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
      do i=1, nlevs
         do k=1, n_points/64
            write(iu_dump, rec=irec) (tac(l, i), l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (tac(l, i)                         &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
!
!     options for treating clouds
      write(iu_dump, rec=irec) global_cloud_top, l_microphysics         &
     &  , l_cloud_water_partition, lcv_3d_cca
      irec=irec+1
!
!     stratiform clouds:
      do i=1, nclds
         do k=1, n_points/64
            write(iu_dump, rec=irec) (lca_area(l, i)                    &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (lca_area(l, i)                    &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
      do i=1, nclds
         do k=1, n_points/64
            write(iu_dump, rec=irec) (lca_bulk(l, i)                    &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (lca_bulk(l, i)                    &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
      do i=1, nclds
         do k=1, n_points/64
            write(iu_dump, rec=irec) (lccwc1(l, i)                      &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (lccwc1(l, i)                      &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
      do i=1, nclds
         do k=1, n_points/64
            write(iu_dump, rec=irec) (lccwc2(l, i)                      &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (lccwc2(l, i)                      &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
!
!     convective clouds:
      do i=1, n_cca_lev
         do k=1, n_points/64
            write(iu_dump, rec=irec) (cca(l, i)                         &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (cca(l, i)                         &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
      do k=1, n_points/64
         write(iu_dump, rec=irec) (cccwp(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (cccwp(l), l=64*(n_points/64)+1       &
     &         , n_points)
         irec=irec+1
      endif
      do k=1, n_points/64
         write(iu_dump, rec=irec) (ccb(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (ccb(l), l=64*(n_points/64)+1         &
     &         , n_points)
         irec=irec+1
      endif
      do k=1, n_points/64
         write(iu_dump, rec=irec) (cct(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (cct(l), l=64*(n_points/64)+1         &
     &         , n_points)
         irec=irec+1
      endif
!
!     surface fields:
      do i=1, 4
        do k=1, n_points/64
           write(iu_dump, rec=irec)                                     &
     &        (land_albedo(l, i), l=64*(k-1)+1, 64*k)
           irec=irec+1
        enddo
        if (n_points >  64*(n_points/64))  then
           write(iu_dump, rec=irec)                                     &
     &        (land_albedo(l, i), l=64*(n_points/64)+1, n_points)
           irec=irec+1
        endif
      enddo
      write(iu_dump, rec=irec) l_moses_ii, l_ctile
      irec=irec+1
      do k=1, n_points/64
         write(iu_dump, rec=irec) (land_alb(l)                          &
     &         , l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (land_alb(l)                          &
     &         , l=64*(n_points/64)+1, n_points)
         irec=irec+1
      endif
      do k=1, n_points/64
         write(iu_dump, rec=irec) (sice_alb(l)                          &
     &         , l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (sice_alb(l)                          &
     &         , l=64*(n_points/64)+1, n_points)
         irec=irec+1
      endif
      do k=1, n_points/64
         write(iu_dump, rec=irec) (flandg(l)                            &
     &         , l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (flandg(l)                            &
     &         , l=64*(n_points/64)+1, n_points)
         irec=irec+1
      endif
      do i=1, 2
         do k=1, n_points/64
            write(iu_dump, rec=irec) (open_sea_albedo(l, i)             &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (open_sea_albedo(l, i)             &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
      do k=1, n_points/64
         write(iu_dump, rec=irec)                                       &
     &      (ice_fraction(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec)                                       &
     &      (ice_fraction(l), l=64*(n_points/64)+1, n_points)
         irec=irec+1
      endif
      do k=1, n_points/64
         write(iu_dump, rec=irec) (land(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (land(l), l=64*(n_points/64)+1        &
     &         , n_points)
         irec=irec+1
      endif
      do k=1, n_points/64
         write(iu_dump, rec=irec) (land0p5(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (land0p5(l), l=64*(n_points/64)+1     &
     &         , n_points)
         irec=irec+1
      endif
      do k=1, n_points/64
         write(iu_dump, rec=irec) (lying_snow(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (lying_snow(l), l=64*(n_points/64)+1  &
     &         , n_points)
         irec=irec+1
      endif
!
!     solar fields:
      do k=1, n_points/64
         write(iu_dump, rec=irec) (coszin(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (coszin(l), l=64*(n_points/64)+1      &
     &         , n_points)
         irec=irec+1
      endif
      do k=1, n_points/64
         write(iu_dump, rec=irec) (lit(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (lit(l), l=64*(n_points/64)+1         &
     &         , n_points)
         irec=irec+1
      endif
      do k=1, nlit/64
         write(iu_dump, rec=irec) (list(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (nlit >  64*(nlit/64))  then
         write(iu_dump, rec=irec) (list(l), l=64*(nlit/64)+1            &
     &         , nlit)
         irec=irec+1
      endif
      write(iu_dump, rec=irec) scs
      irec=irec+1
!
!     aerosols:
      write(iu_dump, rec=irec) l_climat_aerosol, l_clim_aero_hgt        &
     &  , n_levels_bl                                                   &
     &  , l_use_sulpc_direct, l_use_sulpc_indirect                      &
     &  , sulp_dim1, sulp_dim2
      irec=irec+1
      do k=1, n_points/64
         write(iu_dump, rec=irec) (bl_depth(l), l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (bl_depth(l), l=64*(n_points/64)+1    &
     &         , n_points)
         irec=irec+1
      endif
      do i=1, sulp_dim2
         do k=1, sulp_dim1/64
            write(iu_dump, rec=irec) (accum_sulphate(l, i)              &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (sulp_dim1 >  64*(sulp_dim1/64))  then
            write(iu_dump, rec=irec) (accum_sulphate(l, i)              &
     &            , l=64*(sulp_dim1/64)+1, sulp_dim1)
            irec=irec+1
         endif
      enddo
      do i=1, sulp_dim2
         do k=1, sulp_dim1/64
            write(iu_dump, rec=irec) (aitken_sulphate(l, i)             &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (sulp_dim1 >  64*(sulp_dim1/64))  then
            write(iu_dump, rec=irec) (aitken_sulphate(l, i)             &
     &            , l=64*(sulp_dim1/64)+1, sulp_dim1)
            irec=irec+1
         endif
      enddo
      do i=1, sulp_dim2
         do k=1, sulp_dim1/64
            write(iu_dump, rec=irec) (diss_sulphate(l, i)               &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (sulp_dim1 >  64*(sulp_dim1/64))  then
            write(iu_dump, rec=irec) (diss_sulphate(l, i)               &
     &            , l=64*(sulp_dim1/64)+1, sulp_dim1)
            irec=irec+1
         endif
      enddo
!     Seasalt
      write(iu_dump, rec=irec) l_use_seasalt_indirect                   &
     &  , l_use_seasalt_direct, salt_dim_a, salt_dim_b
      irec=irec+1
      do i=1, salt_dim_b
         do k=1, salt_dim_a/64
            write(iu_dump, rec=irec) (sea_salt_film(l, i)               &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (salt_dim_a >  64*(salt_dim_a/64))  then
            write(iu_dump, rec=irec) (sea_salt_film(l, i)               &
     &            , l=64*(salt_dim_a/64)+1, salt_dim_a)
            irec=irec+1
         endif
      enddo
      do i=1, salt_dim_b
         do k=1, salt_dim_a/64
            write(iu_dump, rec=irec) (sea_salt_jet(l, i)                &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (salt_dim_a >  64*(salt_dim_a/64))  then
            write(iu_dump, rec=irec) (sea_salt_jet(l, i)                &
     &            , l=64*(salt_dim_a/64)+1, salt_dim_a)
            irec=irec+1
         endif
      enddo
      write(iu_dump, rec=irec) l_use_soot_direct, soot_dim1, soot_dim2
      irec=irec+1
      do i=1, soot_dim2
         do k=1, soot_dim1/64
            write(iu_dump, rec=irec) (fresh_soot(l, i)                  &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (soot_dim1 >  64*(soot_dim1/64))  then
            write(iu_dump, rec=irec) (fresh_soot(l, i)                  &
     &            , l=64*(soot_dim1/64)+1, soot_dim1)
            irec=irec+1
         endif
      enddo
      do i=1, soot_dim2
         do k=1, soot_dim1/64
            write(iu_dump, rec=irec) (aged_soot(l, i)                   &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (soot_dim1 >  64*(soot_dim1/64))  then
            write(iu_dump, rec=irec) (aged_soot(l, i)                   &
     &            , l=64*(soot_dim1/64)+1, soot_dim1)
            irec=irec+1
         endif
      enddo
!     Murk aerosol
      do i=1, nlevs
         do k=1, n_points/64
            write(iu_dump, rec=irec) (aero_meso(l, i)                   &
     &            , l=64*(k-1)+1, 64*k)
            irec=irec+1
         enddo
         if (n_points >  64*(n_points/64))  then
            write(iu_dump, rec=irec) (aero_meso(l, i)                   &
     &            , l=64*(n_points/64)+1, n_points)
            irec=irec+1
         endif
      enddo
      write(iu_dump, rec=irec) l_murk_rad
      irec=irec+1
      do k=1, n_points/64
         write(iu_dump, rec=irec) (trindx(l)                            &
     &         , l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_points >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (trindx(l)                            &
     &         , l=64*(n_points/64)+1, n_points)
         irec=irec+1
      endif
      write(iu_dump, rec=irec) pts, l_scale_inc                         &
     &  , i_segment, i_call
      irec=irec+1
!
!     Satellite Fields:
      write(iu_dump, rec=irec) n_viewing_direction, n_viewing_level     &
     &  , n_channel
      irec=irec+1
      do i=1, 2
        do j=1, n_viewing_direction
          do k=1, n_points/64
             write(iu_dump, rec=irec) (viewing_direction(l, j, i)       &
     &             , l=64*(k-1)+1, 64*k)
             irec=irec+1
          enddo
          if (n_points >  64*(n_points/64)) then
             write(iu_dump, rec=irec) (viewing_direction(l, j, i)       &
     &             , l=64*(n_points/64)+1, n_points)
             irec=irec+1
          endif
        enddo
      enddo
      do k=1, n_viewing_level/64
         write(iu_dump, rec=irec) (viewing_level(l)                     &
     &         , l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (n_viewing_level >  64*(n_viewing_level/64))  then
         write(iu_dump, rec=irec) (viewing_level(l)                     &
     &         , l=64*(n_viewing_level/64)+1, n_viewing_level)
         irec=irec+1
      endif
      do k=1, nd_band/64
         write(iu_dump, rec=irec) (map_channel(l)                       &
     &         , l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (nd_band >  64*(nd_band/64))  then
         write(iu_dump, rec=irec) (map_channel(l)                       &
     &         , l=64*(nd_band/64)+1, nd_band)
         irec=irec+1
      endif
      do k=1, nlit/64
         write(iu_dump, rec=irec) (row_list(l)                          &
     &         , l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (nlit >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (row_list(l)                          &
     &         , l=64*(nlit/64)+1, nlit)
         irec=irec+1
      endif
      do k=1, nlit/64
         write(iu_dump, rec=irec) (col_list(l)                          &
     &         , l=64*(k-1)+1, 64*k)
         irec=irec+1
      enddo
      if (nlit >  64*(n_points/64))  then
         write(iu_dump, rec=irec) (col_list(l)                          &
     &         , l=64*(nlit/64)+1, nlit)
         irec=irec+1
      endif
      write(iu_dump, rec=irec) l_flux_below_690nm_surf
      irec=irec+1
!
      close(iu_dump)
!
!
!
      return
      END SUBROUTINE r2_sw_dump
#endif
