#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      Subroutine BL_TRACER_INTCTL(                                      &
     &        row_length, rows, bl_levels, off_x, off_y, at_extremity   &
     &,       dtrdz_charney_grid,r_rho_levels,r_theta_levels            &
     &,       halo_i,halo_j                                             &
     &,       tr_levels, tr_vars, model_levels                          &
     &,       DIM_CS2                                                   &
     &,       alpha_cd, rhokh_mix                                       &
     &,       p_star, p, timestep                                       &
! Control logicals
     &,       L_MURK,L_MURK_ADVECT,L_BL_TRACER_MIX,L_DUST,L_CAM_DUST    &
     &,       L_SULPC_SO2, l_sulpc_nh3, l_sulpc_dms, l_soot, l_biomass  &
     &,       l_ocff, l_co2_interactive                                 &
     &,       l_co2_emits                                               &
! rml 2/7/13
     &,       L_CO2_TRACER                                              &
     &,       L_ukca, L_USE_CARIOLLE                                    &
! Fields to mix
     &,       murk, free_tracers, OZONE_TRACER                          &
! Mineral Dust
     &,        DUST_DIV1,DUST_DIV2,DUST_DIV3                            &
     &,        DUST_DIV4,DUST_DIV5,DUST_DIV6                            &
! Sulphur cycle
     &,       so2, dms, so4_aitken, so4_accu, so4_diss, nh3             &
! Soot cycle
     &,       soot_new, soot_aged, soot_cld                             &
! Biomass aerosol
     &,       bmass_new, bmass_agd, bmass_cld                           &
! Fossil-fuel organic carbon aerosol
     &,       ocff_new, ocff_agd, ocff_cld                              &
! Carbon cycle
     &,       co2, co2_emits, co2flux, npp, resp_s                      &
     &,       co2_flux_tot, land_co2                                    &
! Tracer fluxes - kdcorbin, 05/10
     &,       tracer_flux1, tracer_flux2, tracer_flux3, tracer_flux4    &
     &,       tracer_flux5, tracer_flux6, tracer_flux7, tracer_flux8    &
     &,       tracer_flux9, tracer_flux10,tracer_flux11,tracer_flux12   &
     &,       tracer_flux13,tracer_flux14,tracer_flux15,tracer_flux16   &
     &,       tracer_flux17,tracer_flux18,tracer_flux19,tracer_flux20   &     
! Emissions fields
     &, DUST_FLUX                                                       &
     &,       so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em    &
     &,       ocff_hilem, ocff_em, drydep_str                           &
     &,       kent, we_lim, t_frac, zrzi                                &
     &,       kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                &
     &,       zh, zhsc, z_half                                          &
! IN for dry deposition of tracers
     &,       rho_aresist, aresist, resist_b                            &
     &,       R_B_DUST,TSTAR                                            &
     &,       land_points, land_index, ice_frac                         &
! IN variables for MOSES II
     &,       ntype, ntiles, tile_pts, tile_index, tile_frac            &
     &,       canopy, catch, snow_tile, gc                              &
     &,       aresist_tile, resist_b_tile, FLANDG                       &
     &,                                                                 &
! STASH related variables
#include "argsts.h"
     &        stashwork                                                 &
     &        )
!
! ---------------------------------------------------------------------
! Purpose: Intermediate control routine to call BL_TRMIX_DD for tracer
!          variables that require boundary layer mixing and/or dry
!          deposition with MOSES II surface scheme.
!
!          Called by IMP_CTL2.
!
! Current code owner: M Woodage
!
! History:
! Version   Date   Comment
! -------   ----   -------
!  5.3   20/10/01  New deck.                               M Woodage
!  5.4   29/08/02  Extra variables passed through for changes
!                  to tr_mix.                   P. Davison.
!    5.5  26/02/03  Pass 3D CO2 to tracer mixing routine. C Jones
!  5.5   05/02/03  Pass through extra arguments for biomass aerosol
!                  scheme.                           P Davison
!  5.5   19/02/03  Remove redundent RML code. Adrian Lock
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
!   6.2   08/07/05  Added logical and call to copy BL variables to
!                   ukca module.          F. O'Connor
!  6.2  01/03/06  pass CO2 fluxes back up for use in diagnostics
!                 routines.                             C.D. Jones
!  6.1   11/01/06  Remove MOSES I code                A.P.Lock
!
!
!
! Code description:
!   Language: FORTRAN 77  + common extensions
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      Implicit None

#include "c_dust_ndiv.h"
!
! Parallel setup variables
      Integer, Intent(In) ::                                            &
     &  off_x                                                           &
                              ! Size of small halo in i
     &, off_y                                                           &
                              ! Size of small halo in j.
     &, halo_i                                                          &
                              ! Size of halo in i direction
     &, halo_j                ! Size of halo in j direction

       Logical, Intent(In) ::                                           &
     &  at_extremity(4)   ! Indicates if this processor is at north,
                          ! south, east or west of the processor grid

! Model dimensions
      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, bl_levels                                                       &
     &, tr_levels                                                       &
                          ! number of free tracer levels
     &, tr_vars                                                         &
                          ! number of free tracer variables
     &, DIM_CS2         ! soil carbon dimension


! Model switches
      Logical, Intent(In) ::                                            &
     &  L_Murk                                                          &
                          ! Switch for (visibility) aerosol
     &, L_murk_advect                                                   &
                          ! Switch for advecting aerosol
     &, L_bl_tracer_mix                                                 &
                          ! Switch for BL mixing of free tracers
     &, L_DUST                                                          &
                          ! Switch for Mineral Dust
     &, L_CAM_DUST                                                      &
                          !Old version of dust_uplift scheme used in CAM NWP models
     &, L_sulpc_so2                                                     &
                          ! Switch for Sulphur Cycle
     &, L_sulpc_nh3                                                     &
                          ! NH3 included in Sulphur Cycle
     &, L_sulpc_dms                                                     &
                          ! DMS included in Sulphur Cycle
     &, L_soot                                                          &
                          ! Switch for Soot Cycle
     &, L_biomass                                                       &
                          ! Switch for Biomass aerosol
     &, L_ocff                                                          &
                          ! Switch for Fossil-fuel OC aerosol
     &, L_co2_interactive                                               &
                          ! Switch for interactive CO2
     &, L_co2_emits                                                     &
                          ! Switch for anthro CO2 emissions
! rml 2/7/13
     &, L_CO2_tracer                                                    &
                          ! Switch to put co2 flux into passive tracer
     &, L_ukca                                                          &
                          ! switch for UKCA
     &, L_USE_CARIOLLE    
                          ! Switch for the cariolle ozone tracer scheme
!
! Model parameters
      Real, Intent(In) ::                                               &
     &  timestep                                                        &
     &, alpha_cd(bl_levels)                                             &
     &, rhokh_mix (row_length, rows, bl_levels)                         &
     &, drydep_str(row_length, rows)                                    &
     &, we_lim(row_length,rows,3)                                       &
                                      !rho*entrainment rate implied by
                                      !placing of subsidence
     &, zrzi(row_length,rows,3)                                         &
                                      !(z-z_base)/(z_i-z_base)
     &, t_frac(row_length,rows,3)                                       &
                                      !a fraction of the timestep
     &, we_lim_dsc(row_length,rows,3)                                   &
                                      !rho*entrainment rate implied by
                                      !  placing of subsidence
     &, zrzi_dsc(row_length,rows,3)                                     &
                                      !(z-z_base)/(z_i-z_base)
     &, t_frac_dsc(row_length,rows,3)                                   &
                                      !a fraction of the timestep
     &, z_half(row_length,rows,BL_LEVELS)                               &
                                      !Z_HALF(*,K) is height of half
                                      !    level k-1/2.
     &, zhsc(row_length,rows)                                           &
                                      !Top of decoupled layer
     &, zh  (row_length,rows)                                           &
                                      !Top of surface mixed layer
     &, dtrdz_charney_grid(row_length,rows,BL_LEVELS)
                                      ! For tracer mixing

      Integer, Intent(In) ::                                            &
     &  kent(row_length,rows)                                           &
                                      !grid-level of SML inversion
     &, kent_dsc(row_length,rows)                                       &
                                      !grid-level of DSC inversion
     &, land_points                                                     &
                                      !No.of land points being processed
     &, land_index(land_points)                                         &
                                      !set from land_sea_mask
! For MOSES II
     &, ntype                                                           &
                                      !No. of tile types
     &, ntiles                                                          &
                                      !No.of land-surface tiles
     &, tile_pts(ntype)                                                 &
                                      !No.of tile points
     &, tile_index(land_points,ntype) !Index of tile points
!
! Required fields (input)
      Real, Intent(In) ::                                               &
     &  p_star(row_length, rows)                                        &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:bl_levels)              &
                                      ! Heights of theta levels
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, bl_levels)
                                      ! Heights of rho levels


! Mineral Dust Emissions fields (input)
      Real, Intent(In) ::                                               &
     & DUST_FLUX(ROW_LENGTH, ROWS, NDIV)! dust emission flux

! Emissions fields (input)
      Real, Intent(In) ::                                               &
     &  so2_hilem ( row_length, rows )                                  &
     &, so2_em    ( row_length, rows )                                  &
     &, nh3_em    ( row_length, rows )                                  &
     &, dms_em    ( row_length, rows )                                  &
     &, soot_hilem( row_length, rows )                                  &
     &, soot_em   ( row_length, rows )                                  &
     &, ocff_hilem( row_length, rows )                                  &
     &, ocff_em   ( row_length, rows )                                  &
     &, co2_emits ( row_length, rows )                                  &
                                        ! anthro CO2 emissions
     &, co2flux   ( row_length, rows )  ! ocean  CO2 emissions

! Tracer fluxes for atmospheric tracers - kdcorbin, 05/10
      Real, Intent(In) ::                                               &
     &   tracer_flux1(row_length,rows),  tracer_flux2(row_length,rows)   &
     &,  tracer_flux3(row_length,rows),  tracer_flux4(row_length,rows)   &
     &,  tracer_flux5(row_length,rows),  tracer_flux6(row_length,rows)   &
     &,  tracer_flux7(row_length,rows),  tracer_flux8(row_length,rows)   &
     &,  tracer_flux9(row_length,rows),  tracer_flux10(row_length,rows)  &
     &,  tracer_flux11(row_length,rows), tracer_flux12(row_length,rows)  &
     &,  tracer_flux13(row_length,rows), tracer_flux14(row_length,rows)  &
     &,  tracer_flux15(row_length,rows), tracer_flux16(row_length,rows)  &
     &,  tracer_flux17(row_length,rows), tracer_flux18(row_length,rows)  &
     &,  tracer_flux19(row_length,rows), tracer_flux20(row_length,rows)

!
! For dry deposition of tracers (input)
      Real, Intent(In) ::                                               &
     &  rho_aresist(row_length,rows)                                    &
                                       !RHOSTAR*CD_STD*VSHR
     &, aresist(row_length,rows)                                        &
                                       !1/(CD_STD*VSHR)
     &, resist_b(row_length,rows)                                       &
                                       !(1/CH-1/CD_STD)/VSHR
     &, R_B_DUST(ROW_LENGTH,ROWS,NDIV)                                  &
                                       ! surf layer res for dust
     &, TSTAR(ROW_LENGTH,ROWS)                                          &
                               ! temperature
     &, ice_frac(row_length,rows)                                       &
                                       !sea_ice fractn (ice/qrclim.ice)
!
! For MOSES II
     &, tile_frac(land_points,ntiles)                                   &
                                       !fractional coverage for each
!                                       surface tile
     &, canopy(land_points,ntiles)                                      &
                                       !Surface/canopy water (kg/m2)
     &, catch(land_points,ntiles)                                       &
                                       !Surface/canopy water capacity
!                                       of snow-free land tiles (kg/m2)
     &, snow_tile(land_points,ntiles)                                   &
                                       !Snow on tiles (kg/m2).
     &, gc(land_points,ntiles)                                          &
                                       !Stomatal conductance to evapn
!                                       for land tiles (m/s)
     &, aresist_tile(land_points,ntiles)                                &
!                                      ! 1/(CD_STD*VSHR) on land tiles
     &, resist_b_tile(land_points,ntiles)                               &
!                                    !(1/CH-1/CD_STD)/VSHR on land tiles
     &, FLANDG(row_length,rows)                                         &
                                     !land fraction of grid box
!                                       (0 if all sea, 1 if all land)
     &, npp(land_points)                                                &
                                     ! net primary productivity
     &, resp_s(DIM_CS2)                                                 &
                                       ! soil respiration
     &, co2_flux_tot(row_length,rows)                                   &
                                       ! OUT total CO2 flux
     &, land_co2(land_points)          ! OUT terrestrial CO2 flux
! Tracer fields for mixing
      Real, Intent(InOut) ::                                            &
     &  murk        ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,free_tracers( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y,                         &
     &                tr_levels, tr_vars)

      Real, Intent(InOut) ::                                            &
     &  DUST_DIV1   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV2   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV3   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV4   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV5   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV6   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )

      Real, Intent(InOut) ::                                            &
     &  so2         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,dms         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_aitken  ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_accu    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_diss    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,nh3         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )
!
      Real, Intent(InOut) ::                                            &
     &  soot_new    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,soot_aged   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,soot_cld    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_new   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_agd   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_cld   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_new    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_agd    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_cld    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,co2         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,OZONE_TRACER( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          

      Real, Intent(InOut) ::                                            &
     &  stashwork(*)

! Argument comdeck declaration
#include "csubmodl.h"
#include "typsts.h"
!
      External BL_TRMIX_DD

! DEPENDS ON: bl_trmix_dd
      Call BL_TRMIX_DD(                                                 &
     &        row_length, rows, bl_levels, off_x, off_y, at_extremity   &
     &,       dtrdz_charney_grid,r_rho_levels,r_theta_levels            &
     &,       halo_i,halo_j                                             &
     &,       tr_levels, tr_vars, model_levels                          &
     &,       DIM_CS2                                                   &
     &,       alpha_cd, rhokh_mix                                       &
     &,       p_star, p, timestep                                       &
! Control logicals
     &,       L_MURK,L_MURK_ADVECT,L_BL_TRACER_MIX,L_DUST,L_CAM_DUST    &
     &,       L_SULPC_SO2,l_sulpc_nh3, l_sulpc_dms, l_soot, l_biomass   &
     &,       l_ocff, l_co2_interactive                                 &
! rml 2/7/13 added flag for putting co2 flux into passive tracer
     &,       l_co2_emits, L_CO2_TRACER, L_USE_CARIOLLE                 &
! Fields to mix
     &,       murk, free_tracers, OZONE_TRACER                          &
! Mineral Dust
     &,       DUST_DIV1,DUST_DIV2,DUST_DIV3                             &
     &,       DUST_DIV4,DUST_DIV5,DUST_DIV6                             &
! Sulphur cycle
     &,       so2, dms, so4_aitken, so4_accu, so4_diss, nh3             &
! Soot cycle
     &,       soot_new, soot_aged, soot_cld                             &
! Biomass aerosol
     &,       bmass_new, bmass_agd, bmass_cld                           &
! Fossil-fuel organic carbon aerosol
     &,       ocff_new, ocff_agd, ocff_cld                              &
! Carbon cycle
     &,       co2, co2_emits, co2flux, npp, resp_s                      &
     &,       co2_flux_tot, land_co2                                    &
! Tracer fluxes - kdcorbin, 05/10
     &,       tracer_flux1, tracer_flux2, tracer_flux3, tracer_flux4    &
     &,       tracer_flux5, tracer_flux6, tracer_flux7, tracer_flux8    &
     &,       tracer_flux9, tracer_flux10,tracer_flux11,tracer_flux12   &
     &,       tracer_flux13,tracer_flux14,tracer_flux15,tracer_flux16   &
     &,       tracer_flux17,tracer_flux18,tracer_flux19,tracer_flux20   &
! Emissions fields
     &,       DUST_FLUX                                                 &
     &,       so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em    &
     &,       ocff_hilem, ocff_em, drydep_str                           &
     &,       kent, we_lim, t_frac, zrzi                                &
     &,       kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                &
     &,       zh, zhsc, z_half                                          &
! IN for dry deposition of tracers
     &,       rho_aresist, aresist, resist_b                            &
     &,       R_B_DUST,TSTAR                                            &
     &,       land_points, land_index, ice_frac                         &
! IN variables for MOSES II
     &,       ntype, ntiles, tile_pts, tile_index, tile_frac            &
     &,       canopy, catch, snow_tile, gc                              &
     &,       aresist_tile, resist_b_tile, FLANDG                       &
     &,                                                                 &
! STASH related variables
#include "argsts.h"
     &        stashwork                                                 &
     &        )
!
      RETURN
      END SUBROUTINE BL_TRACER_INTCTL
!
#endif
