
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


      Subroutine diagnostics_bl(                                        &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels, boundary_layer_levels    &
     &,                      land_points, dsm_levels                    &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      DIM_CS1, DIM_CS2                           &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      sq_T1p5, L_murk, L_pc2, Aerosol, RHcrit    &
     &,                      L_plsp, plsp, cca_2d, lq_mix_bl            &
     &,                      ls_rain, ls_snow, conv_rain, conv_snow     &
     &,                      at_extremity                               &
     &,                      timestep                                   &
     &,                      p_star                                     &
     &,                      rhcpt, rhc_row_length, rhc_rows            &
     &,                      cloud_fraction_method,overlap_ice_liquid   &
     &,                      ice_fraction_method,ctt_weight,t_weight    &
     &,                      qsat_fixed,sub_cld                         &
     &,                      x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic   &
     &,                      lsp_ei,lsp_fi,lsp_eic,lsp_fic              &
     &,                      ntml, cumulus, L_eacf, u, v, R_u, R_v      &
     &,                      T, rho1, q, qcl, qcf                       &
     &,                      cf, cfl, cff                               &
     &,                      cf_earliest, cfl_earliest, cff_earliest    &
     &,                      T_incr_diagnostic, q_incr_diagnostic       &
     &,                      qcl_incr_diagnostic, qcf_incr_diagnostic   &
     &,                      cf_incr_diagnostic                         &
     &,                      cfl_incr_diagnostic, cff_incr_diagnostic   &
     &,                      u_incr_diagnostic,v_incr_diagnostic        &
     &,                      exner_theta_levels                         &
     &,                      T1p5m, surf_sens_heat_flux, ml_depth       &
     &,                      u10m, v10m, q1p5m                          &
     &,                      e_sea, h_sea, ei                           &
     &,                      sea_ice_htf, sice_mlt_htf                  &
     &,                      snomlt_surf_htf                            &
     &,                      bl_top                                     &
     &,                      bl_type_1,bl_type_2,bl_type_3,bl_type_4    &
     &,                      bl_type_5,bl_type_6,bl_type_7              &
     &,                      fqt, ftl, z0m_gb, z0m_eff, z0h_eff         &
     &,                      rib, surf_latent_heat_flux, taux, tauy     &
     &,                      wind_mixing_energy                         &
     &,                      t_soil, land_index                         &
     &,                      surf_ht_flux                               &
     &,                      surf_ht_flux_land,surf_ht_flux_sice        &
     &,                      rib_ssi,ftl_ssi,e_ssi,ei_sice              &
     &,                      vshr_land,vshr_ssi                         &
     &,                      taux_land,taux_ssi,tauy_land,tauy_ssi      &
     &,                      radnet_sice,flandg                         &
     &,                      ntiles,npft,land_sea_mask,nice             &
     &,                      sil_orog_land,ho2r2_orog,gs,gpp,npp,resp_p &
     &,                      ecan_tile,esoil_tile,gpp_ft,ftl_tile       &
     &,                      npp_ft,resp_p_ft,resp_s,resp_s_tot         &
     &,                      resp_s_tile               & !kdcorbin, 10/10
     &,                      cs,rib_tile,es,ecan,fsmc,radnet_tile       &
     &,                      tstar_tile,canopy,catch,z0m_tile,g_leaf    &
     &,                      t1p5m_tile,q1p5m_tile,le_tile,ei_tile,olr  &
     &,                      epot_tile,tile_frac                        &
     &,                      co2_flux_tot, land_co2                     &
     &,                      l_co2_interactive                          &
     &,                      L_DUST,DUST_FLUX                           &
     &,                      U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE     &
     &,                      DRYDEP2                                    &
! CABLE
     &, TSOIL_TILE,SMCL_TILE,STHF_TILE,SNOW_DEPTH3L,SNOW_MASS3L         &
     &, SNOW_TMP3L,SNOW_RHO3L,SNOW_RHO1L,SNAGE_TILE                     &
! variables required for soil moisture nudging scheme macro
     &,                      rhokh,resfs,chr1p5m,alpha1,ra,wt_ext       &
     &,                      lai_ft,canht_ft,gc                         &
! MGS extra bl vars for UKCA
     &, rhokh_mix, rho_aresist, aresist, resist_b, r_b_dust             &
     &, dtrdz_charney_grid, kent, we_lim, t_frac, zrzi, kent_dsc        &
     &, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &     STASHwork                                                    &
     &  ,model_domain                                                   &
     &  ,sin_theta_longitude, cos_theta_longitude                       &
     &  ,proc_row_group                                                 &
     &  ,BL_diag,l_cable,surf_radflux, CPOOL_TILE, NPOOL_TILE           &
     &  ,PPOOL_TILE, SOIL_ORDER, TRANSP_TILE, GLAI, PHENPHASE)

! Purpose:
!  Calculates diagnostics generated from boundary layer routines
!  (UM section 3).
!
! Method:
!  Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the boundary
! layer routines called previously. Each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing.
!
!  Diagnostics currently available: (in order calculated)
!
! STASH item (section 0)
! ----------------------------------------------------------------------
! 413 sea ice bottom melt by catagories (nice)
! 414 sea ice surface melt by catagories (nice)
!-----------------------------------------------------------------------
!
! STASH item (all section 3 )
! ----------------------------------------------------------------------
! 173 soil respiration (tiles) - kdcorbin, 10/10
! 181 temperature increment across bl + PC2 routines (model levels)
! 182 humidity    increment across bl + PC2 routines (model levels)
! 183 liq con qCL increment across bl + PC2 routines (model levels)
! 184 ice con qCF increment across bl routines (model levels)
! 185 u wind      increment across bl routines (model levels)
! 186 v wind      increment across bl routines (model levels)
! 192 total  cloud fraction increment across bl routines (model levels)
! 193 liquid cloud fraction increment across bl routines (model levels)
! 194 ice    cloud fraction increment across bl routines (model levels)
! 189 liq temp Tl increment across bl routines (model levels)
! 190 tot hum  qT increment across bl routines (model levels)
! ------------------------------------------------------------
! 209 10m u wind (native grid = 'C')
! 225 10m u wind 'B' grid
!    Note: simple linear horizontal interpolation from 'C' to 'B' grid
! 210 10m v wind (native grid = 'C')
! 226 10m v wind 'B' grid
!    Note: simple linear horizontal interpolation from 'C' to 'B' grid
! 227 10m wind speed 'B' grid
! 254 1.5m liquid water temperature = tl calculated in imp_solver
! 255 1.5m total water = qt  (kg water/ kg air) calculated in imp_solver
! 245 1.5m relative humidity [%]
!    Note: (tl,qt at 1.5m) converted to (t,q at 1.5m) using cloud
!          scheme, assuming level 1 rhcrit. No negative rh allowed.
! 236 1.5m temperature
! 237 1.5m specific humidity
! 248 1.5m fog fraction
! 250 1.5m dewpoint
! 251 Silhouette orographic roughness field (A/S)
! 252 Half of peak to trough height of sub-grid orography
! 253 1.5m mist fraction
! 216 heat fluxes on model (rho_) levels
!    Note: o/p 1->boundary_levels but level 1 is not meaningful - set
!          same values as surface field 217.
! 217 surface heat flux
! 234 surface latent heat flux
! 219 surface and model level u wind stress
!    Note: not enabled at vn5.1 since o/p should be on (theta_) levels
!          0->boundary_levels-1, but would be labelled as
!          1->boundary_levels    under current system.
! 220 surface and model level v wind stress
!    Note: as 219.
! 221 magnitude of wind stress on B-grid
! 238 soil temperature on soil levels
!    Note: mdi=-1.e30 for sea points
! 025 mixed layer depth
! 224 wind mixing energy
! 305-310 boundary layer types
! 222 total moisture flux profile on model levels
! 223 total surface moisture flux
! 241 total surface moisture flux per timestep
! 026 effective roughness length for momentum
! 208 RIB - measure of surface stability
! 304 turbulent cloud depth
! 247 visibility
! 232 surface evaporation weighted by leads
! 228 sensible heat flux over open sea
! 229 soil evaporation
! 231 sublimation of sea-ice (accumulation over timestep)
! 201 melting of bottom of sea-ice (GBM)
! 235 melting of top of sea-ice (GBM)
! 256 heat flux through sea ice on catagories
! 257 heat flux melting surface  sea ice on catagories
! 258 heat flux due to melting of snow
!  24 surface temperature
!  49 sea-ice temperature
! 202 heat flux from surface to deep soil level 1
! 296 evap. from soil surface
! 297 evap. from canopy
! 298 sublim. from surface
! 322 TOA outgoing LW radiation
! 259 canopy conductance
! 261 gross primary productivity
! 262 net primary productivity
! 263 plant respiration
! 293 soil respiration
! 287 canopy evap. on tiles
! 288 transpiration + soil evap. on tiles
! 290 sensible heat flux on tiles
! 294 bulk Richardson number on tiles
! 314 net radiation on tiles
! 316 surface temperature on tiles
! 317 tile fractions
! 318 Leaf area indices on vegetated tiles
! 319 Canopy height on vegetated tiles
! 321 canopy water on tiles
! 322 canopy capacity of tiles
! 324 snow adjusted roughness length of tiles
! 328 1.5m temperature over tiles
! 329 1.5m specific humidity over tiles
! 330 latent heat flux on tiles
! 331 sublimation on tiles
! 341 1.5m land mean temperature over tiles
! 342 1.5m land mean specific humidity over tiles
! 50 exchange coefficient for moisture
! 51 combined soil/stomatol/aerodynamic resistance to evaporation
! 52 ratio of coefficients required for 1.5m T calculation
! 53 grad of sat humidity w.r. to temp between surface and level 1
! 54 aerodynamic resistance (s/m)
! 55 cumulative evaporation from soil due to transpiration
! 289 gross primary productivity on PFTs
! 291 net primary productivity on PFTs
! 292 plant respiration on PFTs
! 313 soil moisture availability on PFTs
! 325 leaf turnover rate of PFTs
! 332 heat flux through sea ice on catagories
! 333 heat flux to surface melt of sea ice on catagories
! 464 Obukhov length
! 465 Friction velocity
! 467 Surface buoyancy flux
! 468 Gradient Richardson number
! 466 Convective velocity scale
! 469 Vertical buoyancy gradient
! 470 Modulus of wind shear
! 471 BL Momentum diffusion
! 472 BL heat diffusion
! 473 Turbulent kinetic energy
! 474 x component of orographic stress
! 475 y component of orographic stress

! 334 land potential evaporation rate
! 335 potential evaporation rates on land tiles
! 462 Stomatal conductance
! 476 Combined boundary layer type diagnostic
! ------------------------------------------------------------
! MGS Extra BL variables needed in STASH for UKCA.
! N.B. ZH is already passed in and renamed 'ml_depth'
!  60 rhokh_mix
!  61 RHO_ARESIST (RHOSTAR*CD_STD*VSHR)
!  62 ARESIST [ 1/(CD_STD*VSHR) ]
!  63 RESIST_B (1/CH-1/CD_STD)/VSHR
!  64 DTRDZ_CHARNEY_GRID
!  65 GRID-LEVEL OF SML INVERSION (kent)
!  66 Rho * entrainment rate
!  67 Fraction of the timestep
!  68 zrzi
!  69 GRID-LEVEL OF DSC INVERSION (kent_dsc)
!  70 Rho * entrainment rate (dsc)
!  71 Fraction of the timestep
!  72 zrzi
!  73 ZHSC  Top of decoupled layer
!  74 Surface layer resist for dust div1
!  75 Surface layer resist for dust div2
!  76 Surface layer resist for dust div3
!  77 Surface layer resist for dust div4
!  78 Surface layer resist for dust div5
!  79 Surface layer resist for dust div6
! ------------------------------------------------------------
! 481-484 Individual Pool Soil Respiration

!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Use bl_diags_mod, Only :                                          &
          strnewbldiag

      Implicit None

!C_DUST_NDIV.............................................................
! Description: Contains parameters for mineral dust code
! Current Code Owner: Stephanie Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  5.5      12/02/03  Original Code.   Stephanie Woodward
!
! Declarations:
!
      INTEGER NDIV        ! number of particle size divisions
      PARAMETER (NDIV = 6)
!.....................................................................

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, sq_T1p5                                                         &
     &, L_Murk                                                          &
     &, L_pc2                                                           &
                         ! Use PC2 cloud scheme
     &, L_psd                                                           &
                         ! Use generic ice particle size distribution
     &, L_plsp

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, boundary_layer_levels                                           &
                              ! number of boundary layer levels
     &, land_points                                                     &
     &, dsm_levels                                                      &
     &, ntiles                                                          &
                         ! number of land surface tiles
     &, npft                                                            &
                         ! number of plant funcional types
     &, nice             ! number of seaice catagories

      Logical                                                           &
     & land_sea_mask(row_length,rows) ! land sea mask


      Integer                                                           &
     &  ntml(row_length, rows)   ! Height of diagnosed BL top

      LOGICAL                                                           &
     & lq_mix_bl              ! TRUE if mixing ratios used in
!                             ! boundary layer code

      Logical                                                           &
     &  cumulus(row_length, rows)  ! Logical indicator of convection

      Integer                                                           &
     &  global_row_length                                               &
                            !IN. NUMBER OF points on a global row
     &, global_rows                                                     &
                            !IN. NUMBER OF global rows
     &, me                                                              &
                            !IN. Processor number
     &, halo_i                                                          &
                            !IN. size of large halo in x direction
     &, halo_j                                                          &
                            !IN. size of large halo in y direction
     &, off_x                                                           &
                            !IN. size of small halo in x direction
     &, off_y                                                           &
                            !IN. size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)                                      &
     &, DIM_CS1, DIM_CS2        ! soil carbon dimensions

      Real                                                              &
     & timestep          ! model timestep (seconds)

! Primary Arrays used in all models
      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &      model_levels)                                               &
     &, p_star(row_length,rows)

      Real                                                              &
     &  Aerosol(row_length, rows)                                       &
     &, rhcrit(wet_model_levels)                                        &
     &, T(row_length, rows, model_levels)                               &
     &, rho1(row_length, rows)                                          &
                                       ! Air density at level 1/ kg m-3
     &, q(row_length,rows, wet_model_levels)                            &
     &, qcl(row_length, rows, wet_model_levels)                         &
     &, qcf(row_length, rows, wet_model_levels)                         &
! _star arrays contain theta, q fields before ni_bl_ctl
     &, cf(row_length, rows, wet_model_levels)                          &
     &, cfl(row_length, rows, wet_model_levels)                         &
     &, cff(row_length, rows, wet_model_levels)                         &
! _earliest values contain the values of temperature, water contents
! and cloud fractions before the boundary layer call.
     &, cf_earliest(row_length, rows, wet_model_levels)                 &
     &, cfl_earliest(row_length, rows, wet_model_levels)                &
     &, cff_earliest(row_length, rows, wet_model_levels)                &
! _earliest arrays contain fields at start of imp_ctl
     &, T_incr_diagnostic(row_length, rows, model_levels)               &
     &, q_incr_diagnostic(row_length,rows, wet_model_levels)            &
     &, qcl_incr_diagnostic(row_length, rows, wet_model_levels)         &
     &, qcf_incr_diagnostic(row_length, rows, wet_model_levels)         &
     &, cf_incr_diagnostic(row_length,rows, wet_model_levels)           &
     &, cfl_incr_diagnostic(row_length, rows, wet_model_levels)         &
     &, cff_incr_diagnostic(row_length, rows, wet_model_levels)         &
! _diagnostic contain u,v increments before imp_solver
     &, u_incr_diagnostic(row_length,   rows, model_levels)             &
     &, v_incr_diagnostic(row_length, n_rows, model_levels)             &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, T1p5m(row_length,rows)                                          &
                                 ! IN tl at 1.5m but also
                                 ! OUT workspace as T at 1.5m
     &, u10m(row_length,rows)                                           &
     &, v10m(row_length,n_rows)                                         &
     &, q1p5m(row_length,rows)                                          &
                                 ! IN qt at 1.5m (kg water/ kg air)
                                 ! OUT workspace as q at 1.5m
     &, T1p5m_land(land_points)                                         &
                                 ! Land mean tl at 1.5m but also
     &, q1p5m_land(land_points)                                         &
                                 ! Land mean qt at 1.5m
!                                ! (kg water/ kg air)

     &, surf_sens_heat_flux(row_length,rows)                            &
     &, ftl_ssi(row_length,rows)                                        &
     &, surf_latent_heat_flux(row_length,rows)                          &
     &, ml_depth(row_length,rows)                                       &
     &, taux(row_length,rows,boundary_layer_levels)                     &
     &, tauy(row_length,n_rows,boundary_layer_levels)                   &
     &, taux_land(row_length,rows)                                      &
     &, taux_ssi(row_length,rows)                                       &
     &, tauy_land(row_length,n_rows)                                    &
     &, tauy_ssi(row_length,n_rows)                                     &
     &, radnet_sice(row_length,n_rows)                                  &
     &, flandg(row_length,n_rows)                                       &
     &, wind_mixing_energy(row_length,rows)                             &
     &, t_soil(land_points,dsm_levels)                                  &
     &, surf_ht_flux(row_length, rows)                                  &
     &, surf_ht_flux_land(row_length, rows)                             &
     &, surf_ht_flux_sice(row_length, rows)                             &
!*APL*DIAGS
     &, bl_top(row_length, rows)                                        &
     &, bl_type_1(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if stable
!                                 !     b.l. diagnosed, 0.0 otherwise.
     &, bl_type_2(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if Sc over
!                                 !     stable surface layer diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_3(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if well
!                                 !     mixed b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_4(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer (not over
!                                 !     cumulus) diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_5(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer over cumulus
!                                 !     diagnosed, 0.0 otherwise.
     &, bl_type_6(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if a
!                                 !     cumulus capped b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_7(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if a
!                                 !     shear-dominated b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, ftl(row_length,rows,boundary_layer_levels)                      &
     &, fqt(row_length,rows,boundary_layer_levels)                      &
     &, z0m_gb(row_length, rows)                                        &
     &, z0m_eff(row_length, rows)                                       &
     &, z0h_eff(row_length, rows)                                       &
     &, rib(row_length, rows)                                           &
     &, rib_ssi(row_length, rows)                                       &
     &, e_sea(row_length, rows)                                         &
     &, e_ssi(row_length, rows)                                         &
     &, h_sea(row_length, rows)                                         &
     &, ei(row_length, rows)                                            &
     &, ei_sice(row_length, rows)                                       &
     &, sea_ice_htf(row_length, rows, nice)                             &
     &, sice_mlt_htf(row_length, rows, nice)                            &
     &, snomlt_surf_htf(row_length, rows)                               &
     &, plsp(row_length, rows)                                          &
     &, cca_2d(row_length, rows)                                        &
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)                                     &
     &, vshr_land(row_length, rows)                                     &
     &, vshr_ssi(row_length, rows)
!
! MGS new bl vars for use by UKCA, intent(in):
      Real                                                              &
     &  rhokh_mix (row_length, rows, boundary_layer_levels)             &
     &, rho_aresist(row_length,rows)                                    &
     &, aresist(row_length,rows)                                        &
     &, resist_b(row_length,rows)                                       &
     &, r_b_dust(row_length,rows,ndiv)                                  &
     &, dtrdz_charney_grid(row_length,rows,boundary_layer_levels)       &
     &, we_lim(row_length,rows,3)                                       &
     &, t_frac(row_length,rows,3)                                       &
     &, zrzi(row_length,rows,3)                                         &
     &, we_lim_dsc(row_length,rows,3)                                   &
     &, t_frac_dsc(row_length,rows,3)                                   &
     &, zrzi_dsc(row_length,rows,3)                                     &
     &, zhsc(row_length,rows)
      Integer                                                           &
     &  kent(row_length,rows), kent_dsc(row_length,rows)

! CABLE
      REAL                                                              &
     &  TSOIL_TILE(land_points,NTILES,dsm_levels)                       &
     &, SMCL_TILE(land_points,NTILES,dsm_levels)                        &
     &, STHF_TILE(land_points,NTILES,dsm_levels)                        &
     &, SNOW_DEPTH3L(land_points,NTILES,3)                              &
     &, SNOW_MASS3L(land_points,NTILES,3)                               &
     &, SNOW_TMP3L(land_points,NTILES,3)                                &
     &, SNOW_RHO3L(land_points,NTILES,3)                                &
     &, SNOW_RHO1L(land_points,NTILES)                                  &
     &, SNAGE_TILE(land_points,NTILES)                                  &
     &, TRANSP_TILE(land_points,NTILES)

! CASA CNP
      REAL                                                              &
     &  CPOOL_TILE(land_points,NTILES,10)                               &
     &, NPOOL_TILE(land_points,NTILES,10)                               &
     &, PPOOL_TILE(land_points,NTILES,12)                               &
     &, SOIL_ORDER(land_points)                                         &
     &, GLAI(land_points,NTILES) &
!      INTEGER                                                           &
     &, PHENPHASE(land_points,NTILES)

!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag

      Real                                                              &
     &  sil_orog_land(land_points)                                      &
                                       ! IN silhouette orog roughness
     &, ho2r2_orog(land_points)                                         &
                                       ! IN 1/2 pk to trough ht of orog
     &, gs(land_points)                                                 &
                                       ! IN Canopy conductance
     &, gpp(land_points)                                                &
                                       ! IN Gross primary productivity
     &, npp(land_points)                                                &
                                       ! IN Net primary productivity
     &, resp_p(land_points)                                             &
                                       ! IN Plant respiration
     &, ecan_tile(land_points,ntiles)                                   &
                                       ! IN Canopy evaporation on tiles
     &, esoil_tile(land_points,ntiles)                                  &
                                       ! IN Soil evap. on tiles
     !kdcorbin, 11/10 - changed from npft
     &, gpp_ft(land_points,ntiles)                                      &
                                       ! IN GPP on PFTs
     &, ftl_tile(land_points,ntiles)                                    &
                                       ! IN Sensible heat flux on tiles
     !kdcorbin, 11/10 - changed from npft
     &, npp_ft(land_points,ntiles)                                      &
                                       ! IN NPP on PFTs
     !kdcorbin, 11/10 - changed from npft
     &, resp_p_ft(land_points,ntiles)                                   &
                                       ! IN RESP_P on PFTs
     &, resp_s(land_points,DIM_CS1)                                     &
                                       ! IN Soil respiration
     &, resp_s_tot(DIM_CS2)                                             &
                                       ! IN Soil respiration
     &, resp_s_tile(land_points,ntiles)                                 &
                                       ! IN Soil respiration
                                       ! kdcorbin, 10/10
     &, cs(land_points,DIM_CS1)                                         &
                                       ! IN Soil carbon
     &, cs_tot(DIM_CS2)                                                 &
                                       ! IN diagnosed total soil carbon
     &, rib_tile(land_points,ntiles)                                    &
                                       ! IN RIB on tiles
     &, es(row_length, rows)                                            &
                                       ! IN Evap from soil surface
     &, ecan(row_length, rows)                                          &
                                       ! IN Evap from canopy
     &, fsmc(land_points,npft)                                          &
                                       ! IN Soil moisture availability
!                                      !    factor on PFTs
     &, radnet_tile(land_points,ntiles)                                 &
                                       ! IN Net radiation on tiles
! Lestevens
     &, surf_radflux(row_length,rows)                                   &
                                       ! IN Surface Net Radiation   
     &, tstar_tile(land_points,ntiles)                                  &
                                       ! IN Tile surface temperatures
     &, lai_ft(land_points,npft)                                        &
                                       ! IN LAI of PFTs
     &, canht_ft(land_points,npft)                                      &
                                       ! IN Canopy height of PFTs
     &, gc(land_points,ntiles)                                          &
                                       ! IN Stomatal conductance on tile
     &, canopy(land_points,ntiles)                                      &
                                       ! IN Canopy water on tiles
     &, catch(land_points,ntiles)                                       &
                                       ! IN Canopy capacity on tiles
     &, z0m_tile(land_points,ntiles)                                    &
                                       ! IN Roughness length on tiles
     !kdcorbin, 11/10 - changed from npft
     &, g_leaf(land_points,ntiles)                                      &
                                        ! IN Leaf turnover rate for PFTs
     &, t1p5m_tile(land_points,ntiles)                                  &
                                       ! IN 1.5m temperature over tiles
     &, q1p5m_tile(land_points,ntiles)                                  &
                                       ! IN 1.5m specific humidity
!                                      !    over tiles
     &, le_tile(land_points,ntiles)                                     &
                                       ! IN Surface latent heat flux
!                                      !    on tiles
     &, ei_tile(land_points,ntiles)                                     &
                                       ! IN Sublimation on tiles
     &, olr(row_length, rows)                                           &
                                       ! IN TOA outgoing LW radiation
     &, epot_tile(land_points,ntiles)                                   &
                                       ! IN Potential evap on tiles
     &,tile_frac(land_points,ntiles)                                    &
                              ! IN fractional coverage for each
                              !    surface tile
     &, co2_flux_tot(row_length, rows)                                  &
                                       ! IN total CO2 flux
     &, land_co2(land_points)             ! IN terrestrial CO2 flux

      Logical                                                           &
     & l_co2_interactive                                                &
                                       ! IN Switch for interactive CO2
     &, L_DUST                                                          &
               !IN switch for mineral dust
     &,l_eacf &                         ! IN Empirically adj. cloud frac
     &,l_cable

! Variables for STASH macro for soil moisture nudging scheme
      REAL                                                              &
     & RHOKH(row_length, rows, boundary_layer_levels)                   &
     &,RESFS(land_points, ntiles)                                       &
     &,CHR1P5M(land_points, ntiles)                                     &
     &,ALPHA1(land_points, ntiles)                                      &
     &,RA(land_points)                                                  &
     &,WT_EXT(land_points,dsm_levels)

      Integer                                                           &
     & land_index(land_points)

! Variables required for conversion of 1.5m TL and QT
      Integer                                                           &
     & rhc_row_length, rhc_rows ! rhcpt dimensions
!  rhcpt is point-varying rhcrit and can be a 1-d parameter or a
!  3-d diagnostic field, but only level 1 is required here
      Real                                                              &
     & rhcpt(rhc_row_length, rhc_rows, 1)

! mineral dust variables

      REAL                                                              &
     &  DUST_FLUX(ROW_LENGTH,ROWS,NDIV)                                 &
                                        !dust emission flux
     &, U_S_T_TILE(LAND_POINTS,NTILES,NDIV)                             &
                                           !OUT threshold frict. vel
     &, U_S_T_DRY_TILE(LAND_POINTS,NTILES,NDIV)                         &
                                               !OUT dry soil value
     &, U_S_STD_TILE(LAND_POINTS,NTILES)                                &
                                        !OUT friction velocity
     &, DRYDEP2(ROW_LENGTH,ROWS,NDIV) !dep from 2nd level


! CSUBMODL start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.5    07/04/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!  1. Internal model and submodel dump partition identifiers - fixed
!     for all experiments.
! CSMID start
!
! Description:
!    Hold parameters defining internal model identifiers and submodel
!    data partition (ie main D1 data array and consequent dump), both
!    short and long form.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.3    26/10/93   M. Carter. Part of an extensive mod that:
!                    1.Removes the limit on primary STASH item numbers.
!                    2.Removes the assumption that (section,item)
!                      defines the sub-model.
!                    3.Thus allows for user-prognostics.
!                    Add index to submodel home dump.
! 3.5    13/03/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
! 6.0    02/07/03   Add X_IM and X_SM for small exec.      E.Leung
!
! Declarations:
!
!   Hold parameters defining internal model identifiers and submodel
!   data partition (ie main D1 data array and consequent dump), both
!   short and long form
      ! Internal models
      INTEGER,PARAMETER:: A_IM      = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: ATMOS_IM  = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: O_IM      = 2 ! Ocean internal model
      INTEGER,PARAMETER:: OCEAN_IM  = 2 ! Ocean internalmodel
      INTEGER,PARAMETER:: S_IM      = 3 ! Slab internal model
      INTEGER,PARAMETER:: SLAB_IM   = 3 ! Slab internal model
      INTEGER,PARAMETER:: W_IM      = 4 ! Wave internal model
      INTEGER,PARAMETER:: WAVE_IM   = 4 ! Wave internal model
      INTEGER,PARAMETER:: I_IM      = 5 ! Sea=ice internal model
      INTEGER,PARAMETER:: SEAICE_IM = 5 ! Sea=ice internal model
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_IM      = 6 ! ND internal model
      INTEGER,PARAMETER:: NATMOS_IM = 6 ! ND internal model
      ! Small Executables
      INTEGER,PARAMETER:: X_IM      = 7 ! SX indicator

      ! Submodels
      INTEGER,PARAMETER:: A_SM      = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: ATMOS_SM  = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: O_SM      = 2 ! Ocean submodel
      INTEGER,PARAMETER:: OCEAN_SM  = 2 ! Ocean submodel
      INTEGER,PARAMETER:: W_SM      = 4 ! Wave submodel
      INTEGER,PARAMETER:: WAVE_SM   = 4 ! Wave submodel
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_SM      = 6 ! ND submodel
      INTEGER,PARAMETER:: NATMOS_SM = 6 ! ND submodel
      ! Small Executables
      INTEGER,PARAMETER:: X_SM      = 7 ! SX indicator

! CSMID end

!
!  2. Maximum internal model/submodel array sizes for this version.
!
! CSUBMAX start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.5    13/07/95   Original code. D.M. Goddard
! 4.0     3/11/95   Reduce max internal model, submodel from 10 to 4
!                   to save space in model. At 4.0 the max no of
!                   supported models is 3, 1 slot is reserved for
!                   expansion. Rick Rawlins.
!  4.1  21/02/96  Wave model introduced as 4th sub-model.  RTHBarnes
!
! Declarations:
!
!
!  1. Maximum internal model/submodel array sizes for this version.
!
      ! Max no. of internal models
      INTEGER,PARAMETER:: N_INTERNAL_MODEL_MAX=4

      ! Max no. of submodel dump partitions
      INTEGER,PARAMETER:: N_SUBMODEL_PARTITION_MAX=4

      ! Max value of internal model id
      INTEGER,PARAMETER:: INTERNAL_ID_MAX=N_INTERNAL_MODEL_MAX

      ! Max value of submodel dump id
      INTEGER,PARAMETER:: SUBMODEL_ID_MAX=N_SUBMODEL_PARTITION_MAX

! CSUBMAX end
!
!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.
      INTEGER :: N_INTERNAL_MODEL          ! No. of internal models
      INTEGER :: N_SUBMODEL_PARTITION      ! No. of submodel partitions

      ! Internal models
      INTEGER :: INTERNAL_MODEL_LIST(N_INTERNAL_MODEL_MAX)

      ! Submodel identifier for each internal model in list
      INTEGER :: SUBMODEL_FOR_IM    (N_INTERNAL_MODEL_MAX)

      ! Submodel number for each submodel id
      INTEGER :: SUBMODEL_FOR_SM(N_INTERNAL_MODEL_MAX)

      ! Namelist for information in 3.
      NAMELIST/NSUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,          &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM

      ! 4. Lists calculated in model from user interface supplied arrays
      ! experiment specific.

      ! No of internal models in each submodel partition indexed by sm
      !  identifier
      INTEGER :: N_INTERNAL_FOR_SM(SUBMODEL_ID_MAX)

      ! List of  submodel partition identifiers
      INTEGER :: SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION_MAX)

      ! Submodel partition identifier indexed by internal model identifie
      INTEGER :: SUBMODEL_PARTITION_INDEX(INTERNAL_ID_MAX)

      ! Sequence number of internal model indexed by internal model
      ! identifier: required to map from id to STASH internal model
      ! sequence
      INTEGER :: INTERNAL_MODEL_INDEX(INTERNAL_ID_MAX)


      ! Last internal model within a submodel partition if .TRUE.,
      ! indexed by internal model id.
      LOGICAL :: LAST_IM_IN_SM(INTERNAL_ID_MAX)

      ! Common block for information in 3. and 4.
      COMMON/SUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,             &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM,SUBMODEL_FOR_SM,            &
     &  N_INTERNAL_FOR_SM,SUBMODEL_PARTITION_LIST,                      &
     &  SUBMODEL_PARTITION_INDEX,                                       &
     &  INTERNAL_MODEL_INDEX,                                           &
     &  LAST_IM_IN_SM

!
!  5. Time information specifying coupling frequencies between internal
!     models and submodels, and multipliers, indexed by sequence of
!     internal models and submodels (ie left to right along node tree).
!     {Not required at this release}.
!
! Namelists for information in 5. {Not required at this release}
!
!
!  6. Lists of coupling nodes defining coupling frequencies between
!     internal models and between submodel partitions. (Not defined
!     yet at this release).
!CALL CNODE
!
!  7. Variables dealing with general coupling switches at the control
!     level. {These will require revision at the next release when
!     coupling between internal models is dealt with more generally.
!     Logicals below are set in routine SETGRCTL.}

      ! new internal model next group of timesteps if .true.
      LOGICAL :: new_im

      ! new submodel dump  next group of timesteps if .true.
      LOGICAL :: new_sm

      COMMON/CSUBMGRP/new_im,new_sm

      INTEGER SUBMODEL_IDENT
      COMMON/SUBMODID/SUBMODEL_IDENT
! CSUBMODL end
! TYPSTS starts
! CSUBMODL must be included before this file
!Applicable to all configurations (except MOS variables)
!STASH related variables for describing output requests and space
!management.
!LL
!LL   AUTHOR            Rick Rawlins
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL   3.2             Code creation for Dynamic allocation
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL   3.5  Apr. 95   Sub-Models project.
!LL                  Dimensioning of various STASH arrays altered in
!LL                  accordance with internal model separation scheme.
!LL                  Arrays PPXREF, INDEX_PPXREF deleted as they are no
!LL                  longer required.
!LL                  S.J.Swarbrick
!LL
!
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS

      INTEGER :: MOS_MASK_LEN         ! Size of bit mask for MOS

      COMMON/DSIZE_AO/  MOS_MASK_LEN

! TYPSTSZ end
!LL  Comdeck: CPPXREF --------------------------------------------------
!LL
!LL  Purpose: Holds PARAMETER definitions to describe the structure of
!LL           each STASHmaster file record plus some valid entries.
!LL
!LL  Author    Dr T Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                  Add a PPXREF record for model number.
!LL  4.0   26/07/95  T.Johns.  Add codes for real/int/log data types.
!LL  3.5   10/3/94   Sub-Models project:
!LL                 List of PPXREF addressing codes augmented, in order
!LL                 to include all of the pre_STASH master information
!LL                 in the new PPXREF file.
!LL                 PPXREF_CODELEN increased to 38.
!LL                 PPXREF_IDLEN deleted - no longer relevant.
!LL                   S.J.Swarbrick
!LL  4.1   June 96  Wave model parameters included.
!LL                 ppx_ address parameters adjusted to allow for
!LL                  reading option code as 4x5 digit groups.
!LL                   S.J.Swarbrick
!LL  5.0   29/06/99  Add halo type parameter for new dynamics.
!LL                  New grid codes for LAM boundary conditions
!LL                  D.M. Goddard
!LL  5.1   07/03/00  Fixed/Free format conversion
!LL  5.2   19/09/00  Added ppx_atm_lbc_orog descriptor   P.Burton
!LL  5.3   21/08/01  Added ocean lbc descriptors.   M. J. Bell
!LL  5.3   23/07/01  Add valid pp_lbvc codes referenced in UM. R Rawlins
!LL  5.5   30/01/03  Option code increase from 20 to 30 digits thus
!LL                  requiring option code address range increase by
!LL                  2 so all subsequent addressing codes need to be
!LL                  increased by 2 to make a gap.
!LL                  W Roseblade
!LL
!LL  Logical components covered: C40
!LL
!-----------------------------------------------------------------------
! Primary file record definition
      ! length of ID in a record
      Integer, Parameter :: PPXREF_IDLEN      = 2

      ! total length of characters *WARNING* must be multiple of 4
      ! to avoid overwriting
      Integer, Parameter :: PPXREF_CHARLEN    = 36

      ! number of packing profiles
      Integer, Parameter :: PPXREF_PACK_PROFS = 10

      ! total length of codes = no. of codes (excluding profs)
      ! + pack_profs
      Integer, Parameter :: PPXREF_CODELEN    = 33 + PPXREF_PACK_PROFS

! Derived file record sizes
      ! Assume that an integer is at least 4 bytes long. Wastes some
      ! space on an 8 byte machine.
      ! ppx_charword = 9.
      Integer, Parameter :: PPX_CHARWORD      = ((PPXREF_CHARLEN+3)/4)

      ! read buffer record length
      Integer, Parameter :: PPX_RECORDLEN = PPX_CHARWORD+PPXREF_CODELEN
!
!-----------------------------------------------------------------------
! Addressing codes within PPXREF
      Integer, Parameter ::  ppx_model_number   = 1  ! Model number
                                                     ! address
      Integer, Parameter ::  ppx_section_number = 2  ! Section number
                                                     ! address
      Integer, Parameter ::  ppx_item_number    = 3  ! Item number
                                                     ! address
      Integer, Parameter ::  ppx_version_mask   = 4  ! Version mask
                                                     ! address
      Integer, Parameter ::  ppx_space_code     = 5  ! Space code
                                                     ! address
      Integer, Parameter ::  ppx_timavail_code  = 6  ! Time availability
                                                     !  code  address
      Integer, Parameter ::  ppx_grid_type      = 7  ! Grid type code
                                                     ! address
      Integer, Parameter ::  ppx_lv_code        = 8  ! Level type code
                                                     ! address
      Integer, Parameter ::  ppx_lb_code        = 9  ! First level code
                                                     !  address
      Integer, Parameter ::  ppx_lt_code        =10  ! Last level code
                                                     ! address
      Integer, Parameter ::  ppx_lev_flag       =11  ! Level compression
                                                     !  flag  address
      Integer, Parameter ::  ppx_opt_code       =12  ! Sectional option
                                                     ! code  address
      Integer, Parameter ::  ppx_pt_code        =18  ! Pseudo dimension
                                                     ! type  address
      Integer, Parameter ::  ppx_pf_code        =19  ! First pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_pl_code        =20  ! Last pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_ptr_code       =21  ! Section 0 point-
                                                     ! back code address
      Integer, Parameter ::  ppx_dump_packing   =22  ! Dump packing code
                                                     ! address
      Integer, Parameter ::  ppx_lbvc_code      =23  ! PP LBVC code
                                                     ! address
      Integer, Parameter ::  ppx_rotate_code    =24  ! Rotation code
                                                     ! address
      Integer, Parameter ::  ppx_field_code     =25  ! PP field code
                                                     ! address
      Integer, Parameter ::  ppx_user_code      =26  ! User code address
      Integer, Parameter ::  ppx_meto8_levelcode=27  ! CF level code
                                                     ! address
      Integer, Parameter ::  ppx_meto8_fieldcode=28  ! CF field code
                                                     ! address
      Integer, Parameter ::  ppx_cf_levelcode   =27
      Integer, Parameter ::  ppx_cf_fieldcode   =28
      Integer, Parameter ::  ppx_base_level     =29  ! Base level code
                                                     ! address
      Integer, Parameter ::  ppx_top_level      =30  ! Top level code
                                                     ! address
      Integer, Parameter ::  ppx_ref_lbvc_code  =31  ! Ref level LBVC
                                                     ! code address
      Integer, Parameter ::  ppx_data_type      =32  ! Data type code
                                                     ! address
      Integer, Parameter ::  ppx_halo_type      =33
      Integer, Parameter ::  ppx_packing_acc    =34  ! Packing accuracy
                                                     ! code  address
      Integer, Parameter ::  ppx_pack_acc       =34  ! Must be last:


                                                 ! multiple pack_acc to
                                                 ! fill up remaining
                                                 ! array elements


!-------------------------------------------------------------------
! Valid grid type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_atm_nonstd=0      ! Non-standard atmos
                                                  ! grid
      Integer, Parameter :: ppx_atm_tall=1        ! All T points (atmos)
      Integer, Parameter :: ppx_atm_tland=2       ! Land-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tsea=3        ! Sea-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tzonal=4      ! Zonal field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_tmerid=5      ! Merid field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_uall=11       ! All u points (atmos)
      Integer, Parameter :: ppx_atm_uland=12      ! Land-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_usea=13       ! Sea-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_uzonal=14     ! Zonal field at u
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_umerid=15     ! Merid field at u
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_scalar=17     ! Scalar (atmos)
      Integer, Parameter :: ppx_atm_cuall=18      ! All C-grid (u)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_cvall=19      ! All C-grid (v)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_compressed=21 ! Compressed land
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_ozone=22      ! Field on ozone
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_river=23      ! River routing
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_rim=25        ! Rim type field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_theta=26  ! All T points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_u=27      ! All u points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_v=28      ! All v points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_orog=29   ! Orography field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_ocn_nonstd=30     ! Non-standard ocean
                                                  ! grid
      Integer, Parameter :: ppx_ocn_tcomp=31      ! Compressed T points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_ucomp=32      ! Compressed u points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_tall=36       ! All T points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_uall=37       ! All u points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_cuall=38      ! All C-grid (u)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_cvall=39      ! All C-grid (v)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_tfield=41     ! All non-cyclic T
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_ufield=42     ! All non-cyclic u
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_tzonal=43     ! Zonal n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_uzonal=44     ! Zonal n-c field at
                                                  ! u points (ocean)
      Integer, Parameter :: ppx_ocn_tmerid=45     ! Merid n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_umerid=46     ! Merid n-c field at
                                                  ! u points  (ocean)
      Integer, Parameter :: ppx_ocn_scalar=47     ! Scalar (ocean)
      Integer, Parameter :: ppx_ocn_rim=51        ! Rim type field
                                                  ! (LAM BCs ocean)
      Integer, Parameter :: ppx_ocn_lbc_theta=52  ! Ocean rim fields
      Integer, Parameter :: ppx_ocn_lbc_u=53      ! on T & U grids
      Integer, Parameter :: ppx_wam_all=60        ! All points (wave
                                                  ! model)
      Integer, Parameter :: ppx_wam_sea=62        ! Sea points only
                                                  ! (wave model)
      Integer, Parameter :: ppx_wam_rim=65        ! Rim type field
                                                  ! (LAM BCs wave)

!--------------------------------------------------------------------
! Valid rotation type codes
!--------------------------------------------------------------------
      Integer, Parameter :: ppx_unrotated=0       ! Unrotated output
                                                  ! field
      Integer, Parameter :: ppx_elf_rotated=1     ! Rotated ELF field

!-------------------------------------------------------------------
! Valid level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_full_level=1      ! Model full level
      Integer, Parameter :: ppx_half_level=2      ! Model half level
      Integer, Parameter :: ppx_rho_level=1       ! Model rho level
      Integer, Parameter :: ppx_theta_level=2     ! Model theta level
      Integer, Parameter :: ppx_single_level=5    ! Model single level
      Integer, Parameter :: ppx_soil_level=6      ! Deep Soil level

!-------------------------------------------------------------------
! Valid data type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_type_real=1       ! Real data type
      Integer, Parameter :: ppx_type_int=2        ! Integer data type
      Integer, Parameter :: ppx_type_log=3        ! Logical data type

!-------------------------------------------------------------------
! Valid meto8 level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_meto8_surf=9999   ! MetO8 surface type
                                                  ! code

!-------------------------------------------------------------------
! Valid dump packing codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_pack_off=0        ! Field not packed
                                                  ! (ie. 64 bit)
      Integer, Parameter :: ppx_pack_32=-1        ! Field packed to
                                                  ! 32 bit in  dump
      Integer, Parameter :: ppx_pack_wgdos=1      ! Field packed by
                                                  ! WGDOS method
      Integer, Parameter :: ppx_pack_cfi1=11      ! Field packed using
                                                  ! CFI1  (ocean)

!-------------------------------------------------------------------
! Add valid lbvc codes referenced in model (pp header output labels)
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_lbvc_height  =  1 ! height
      Integer, Parameter :: ppx_lbvc_depth   =  2 ! depth (ocean)
      Integer, Parameter :: ppx_lbvc_pressure=  8 ! pressure
      Integer, Parameter :: ppx_lbvc_theta   = 19 ! potential T
      Integer, Parameter :: ppx_lbvc_hybrid  = 65 ! hybrid height(atmos)
      Integer, Parameter :: ppx_lbvc_PV      = 82 ! potential vorticity
      Integer, Parameter :: ppx_lbvc_surface =129 ! surface
! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      INTEGER :: MOS_OUTPUT_LENGTH
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP,MOS_OUTPUT_LENGTH

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
      INTEGER:: MOS_MASK(MOS_MASK_LEN)
! TYPSTS end
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
!  Global Variables:----------------------------------------------------
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
! C_VISBTY start
!LL Description:
!LL   This COMDECK contains declarations for constants used to diagnose
!LL visibility. Constants are set as PARAMTERs.
!LL
!LL
!LL  Model            Modification history:
!LL version  Date
!LL  3.2    29/04/93  CCN Parameters moved here from VISBTY so that
!LL                   they can also be used to compute fog fraction.
!LL                   Programmer: Pete Clark.
!LL  4.0 05/09/95  Variable AEROMAX used as upper limit to aerosol in
!LL                assimilation introduced. Programmer Pete Clark.
!LL  4.5 01/05/98  Completely re-written for NIMROD style diagnostic.
!    5.3 27/07/01  New parameters for the effect of precipitation
!                  visibility created.                 Pete Clark
!LL  5.3 17/10/01  Rename rho and rho_a. Adam Clayton
!LL
!LLEND----------------------------------------------------------------
      INTEGER, PARAMETER :: n_vis_thresh = 2

      ! Standard number density of the aerosol (/m3)
      REAL, PARAMETER :: N0 = 500.0E6

      ! Activation parameter
      REAL, PARAMETER :: B0= 0.5

      ! Radius of standard aerosol particle (m)
      REAL, PARAMETER :: radius0 = 0.16E-6

      REAL, PARAMETER :: FourThirds = 4.0/3.0     ! 4/3

      REAL :: vis_thresh(n_vis_thresh)
      DATA vis_thresh /1000.0,5000.0/
      ! Density of the the aerosol (Kg/m3)
      REAL, PARAMETER :: rho_aerosol = 1700.0

      ! Density of air (Kg/m3)
      REAL, PARAMETER :: rho_air = 1.0
      ! Standard aerosol mass mixing ratio (Kg/Kg)
      REAL, PARAMETER :: m0 = FourThirds * Pi *                         &
     &  radius0 * radius0 * radius0 *                                   &
     &                         (rho_aerosol/rho_air) * N0

      ! Aerosol particle radius/mass loading power
      REAL, PARAMETER :: power = 1.0/6.0

      ! Scattering coefficient normalisation
      REAL, PARAMETER :: Beta0  = 1.5 * Pi

      REAL, PARAMETER :: LiminalContrast   = 0.02

      ! Natural log of Liminal contrast
      REAL, PARAMETER :: LnLiminalContrast = -3.912023005

      ! Constant incorporating the scattering coefficient, normalisation
      ! transformation to visibility ( = ln(liminal contrast) / Beta0 )
      REAL, PARAMETER :: VisFactor = -LnLiminalContrast / Beta0

      ! Reciprocal of the clean air visibility
      REAL, PARAMETER :: RecipVisAir = 1.0E-5

      ! Constant involving surface energy of water
      REAL, PARAMETER :: A0 = 1.2E-9

      ! Visibility defining fog
      REAL, PARAMETER :: VISFOG = 1000.0

      ! Visibility defining mist
      REAL, PARAMETER :: VISMIST = 5000.0

      ! Minimum allowed aerosol
      REAL, PARAMETER :: AERO0 = 0.1

      ! maximum allowed aerosol
      REAL, PARAMETER :: AEROMAX = 200.0

      ! tunable parameter: the cumulative prob value at which vis is
      ! estimated
      REAL, PARAMETER :: calc_prob_of_vis = 0.4
!*L------------------COMDECK CCARBON------------------------------------
! Purpose: declares variables and parameters for the carbon cycle
! History:
! version  date         change
! 5.5      26/02/03     add M_CARBON. C Jones.
!----------------------------------------------------------------------
!carbon cycle and vegetation parameters
      REAL                                                              &
     & M_CO2                                                            &
                                  ! molecular weight of CO2
     &,M_AIR                                                            &
                                  ! molecular weight of dry air
     &,M_CARBON                                                         &
                                  ! molecular weight of carbon
     &,EPSILON                                                          &
                                  ! Ratio of molecular weights of water
!                                 !  and dry air.
     &,EPCO2                                                            &
                                  ! Ratio of molecular weights of CO2
!                                 !  and dry air.
     &,EPO2                                                             &
                                  ! Ratio of molecular weights of O2
!                                 !  and dry air.
     &,CO2CONV_A2O                                                      &
                                  ! conversion factor for atmos to
!                                 !  ocean passing of CO2 (mmr to ppmv)
     &,CO2CONV_O2A                ! conversion factor for ocean to
!                                 !  atmos passing of CO2 flux
!                                 !  (mol C/m2/yr to Kg CO2/m2/s)

      PARAMETER (M_AIR=28.966, EPCO2=1.5194, M_CO2=M_AIR*EPCO2,         &
     &           M_CARBON = 12.0, EPSILON = 0.62198, EPO2 = 1.106)

      PARAMETER (CO2CONV_A2O = M_AIR * 1E6 / M_CO2,                     &
     &           CO2CONV_O2A = M_CO2 * 1e-3 / (360.0 * 24.0 * 3600.0))
!*----------------------------------------------------------------------
!
!  Local parameters and other physical constants------------------------
      REAL LCRCP                        ! Derived parameter.
      PARAMETER ( LCRCP=LC/CP )         ! Lat ht of condensation/Cp.
      Logical, PARAMETER :: l_mixing_ratio=.false. ! Use mixing ratios
!
! Diagnostics info
      Real                                                              &
     & STASHwork(*)                                                     &
                        ! STASH workspace

     &,epot_land(land_points) ! Land mean Potential evap


! Local variables

      Integer                                                           &
     &  i, j, k, l, n                                                   &
                         ! loop indices
     &,    icode                                                        &
                                ! Return code  =0 Normal exit  >1 Error
     & ,item                                                            &
                        ! STASH item
     & ,sect            ! STASH section
      Integer vcode
      Parameter( sect = 3 ) ! for boundary layer

      Character*80  cmessage

      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_bl')

      Integer                                                           &
     &  im_index        ! internal model index

     Integer :: land_points_dum
!                Dummy number of land points returned from the
!                routine FROM_LAND_POINTS: the variable land_points
!                is already set on input

      Integer                                                           &
     &  PP_code_bl_types(7)

      Integer                                                           &
     & PSLEVEL                                                          &
                     !  loop counter for pseudolevels
     &,PSLEVEL_OUT   !  index for pseudolevels sent to STASH

      Logical                                                           &
     & PLLTILE(NTILES)                                                  &
                          ! pseudolevel list for surface types
     &,PLLPFT(NPFT)                                                     &
                          ! pseudolevel list for PFTs
     &,PLLNICE(NICE)                                                    &
     &,PLLBL(npft)
                    ! pseudolevel list for BL diagnostics with 3rd dim=3

      Real                                                              &
     &  interp_data(row_length,rows)                                    &
     &, interp_data_3(row_length*rows*model_levels)                     &
                                                    ! work array
     &, interp_data_bl(row_length,rows,npft)                            &
! work array, for diagnostics with 3rd dim = 3. Arbitrarily use
!  npft, which equals 5.
     &, qcl1p5m(row_length,rows)                                        &
     &, Beta_LS_Rain(row_length, rows)                                  &
                                       ! Scattering in LS Rain.
     &, Beta_LS_Snow(row_length, rows)                                  &
                                       ! Scattering in LS Snow.
     &, Beta_C_Rain(row_length, rows)                                   &
                                       ! Scattering in Conv Rain
     &, Beta_C_Snow(row_length, rows)                                   &
                                       ! Scattering in Conv Snow
     &, Vis(row_length, rows)                                           &
     &, Vis_no_precip(row_length, rows)                                 &
     &, Vis_LS_Precip(row_length, rows)                                 &
     &, Vis_C_Precip(row_length, rows)                                  &
     &, vis_pr(row_length, rows,1,n_vis_thresh)                         &
     &, Vis_Threshold(row_length,rows,1,n_vis_thresh)                   &
                                                     !FOG_FR works for n
                                                     ! levels, we want 1
     &, PVis(row_length,rows,n_vis_thresh)                              &
     &, u10mB(row_length,n_rows)                                        &
                                   !10m  u-wind     B-grid
     &, v10mB(row_length,n_rows)                                        &
                                   !10m  v-wind     B-grid
     &, ws10mB(row_length,n_rows)                                       &
                                   !10m  wind speed B-grid
     &, ws10m_p(row_length,n_rows)                                      &
                                   !10m  wind speed P-grid
     &, tauxB(row_length,n_rows,boundary_layer_levels)                  &
                                                      !x-stress   B-grid
     &, tauyB(row_length,n_rows,boundary_layer_levels)                  &
                                                      !y-stress   B-grid
     &, tauB(row_length,n_rows,boundary_layer_levels)                   &
                                                      !stress mag B-grid
     &, SEA_ICE_HTF_GB(ROW_LENGTH,ROWS)                                 &
                                         ! Gridbox mean flux through ice
     &, SICE_MLT_HTF_GB(ROW_LENGTH,ROWS)                                &
                                         ! Gridbox_mean ice surface melt
     &, bl_type_comb(row_length,rows) ! Combined boundary layer
!                                     ! type diagnostic.

! Variables required for calculating the maximum wind gust
! diagnostic
      real, allocatable :: gust_wind(:,:)
      real, allocatable :: std_dev(:,:)
      real, allocatable :: u10m_halo(:,:)
      real, allocatable :: u10m_p(:,:)
      real, allocatable :: v10m_halo(:,:)
      real, allocatable :: v10m_p(:,:)


! The following are declared as arrays rather than scalars to maintain
! consistency of argument type, since they are passed to a general 
! purpose subroutine in which they are array arguments
      real :: mag_vector_np(1)
      real :: mag_vector_sp(1)
      real :: dir_vector_np(1)
      real :: dir_vector_sp(1)


! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
      real :: scale
      real :: cos_theta_longitude(row_length,rows)
      real :: sin_theta_longitude(row_length,rows)

! Tunable parameters used in the calculation of the maximum
! wind gust. See http://www-nwp/~frpz/gust_diag/max_gust.html
      real, parameter :: c_ugn=4.0
      real, parameter :: gust_const=2.29

      integer :: i_g,j_g,model_domain
      integer :: y,z
      integer :: proc_row_group
      integer, parameter :: PNorth=1
      integer, parameter :: PEast =2
      integer, parameter :: PSouth=3
      integer, parameter :: PWest =4
      integer, parameter :: NoDomain = -1
! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end


      Real                                                              &
                  !, INTENT(IN)
     &  x1i                                                             &
                  ! Intercept of aggregate size distribution
     &, x1ic                                                            &
                  ! Intercept of crystal size distribution
     &, x1r                                                             &
                  ! Intercept of raindrop size distribution
     &, x2r                                                             &
                  ! Scaling parameter of raindrop size distribution
     &, x4r                                                             &
                  ! Shape parameter of raindrop size distribution
     &, ai, bi, aic, bic                                                &
                  ! Mass diameter relationships for ice m(D)=ai D^bi
     &, lsp_ei, lsp_fi, lsp_eic, lsp_fic
                  ! Ice particle Best-Reynolds number relationships
                  ! Re(D) =LSP_EI(C) Be^LSP_FI(C)
      Integer                                                           &
                               !, INTENT(IN)
     &  cloud_fraction_method                                           &
                               ! Method for calculating
                               ! total cloud fraction
     &, ice_fraction_method    ! Method for calculating ice cloud frac.

      Real                                                              &
                               !, INTENT(IN)
     &  overlap_ice_liquid                                              &
                               ! Overlap between ice and liquid phases
     &, ctt_weight                                                      &
                               ! Weighting of cloud top temperature
     &, t_weight                                                        &
                               ! Weighting of local temperature
     &, qsat_fixed                                                      &
                               ! Fixed value of saturation humidity
     &, sub_cld                ! Scaling parameter

! External routines
      External                                                          &
     &  copydiag, copydiag_3d                                           &
     &  ,Ereport                                                        &
     & ,beta_precip,calc_vis_prob,visbty,vis_precip


      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! initialize bl_types pp codes :
      PP_code_bl_types(1)      = 305
      PP_code_bl_types(2)      = 306
      PP_code_bl_types(3)      = 307
      PP_code_bl_types(4)      = 308
      PP_code_bl_types(5)      = 309
      PP_code_bl_types(6)      = 310
      PP_code_bl_types(7)      = 340


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output :
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

!   Copy diagnostic information to STASHwork for STASH processing
!
! ----------------------------------------------------------------------
! DIAG.03181 Copy T   INC: bdy layer + PC2 condensation to stashwork
! ----------------------------------------------------------------------
      item = 181
! Diag03181_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Diag03181_do1:
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              l = i + (j-1)*row_length + (k-1)*rows*row_length
              interp_data_3(l) = T_incr_diagnostic(i, j, k)
            End Do
          End Do
        End Do  ! Diag03181_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       interp_data_3,                                             &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
        cmessage="diag_bl  : error in copydiag_3d(T   INC: bdy lay+PC2)"
        End if
      End if  ! Diag03181_if1
!
! ----------------------------------------------------------------------
! DIAG.03182 Copy q   INC: bdy layer + PC2 condensation to stashwork
! ----------------------------------------------------------------------
      item = 182
! Diag03182_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Diag03182_do1:
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              l = i + (j-1)*row_length + (k-1)*rows*row_length
              interp_data_3(l)  =   q_incr_diagnostic(i, j, k)
            End Do
          End Do
        End Do  ! Diag03182_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       interp_data_3,                                             &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
        cmessage="diag_bl  : error in copydiag_3d(q   INC: bdy lay+PC2)"
        End if
      End if  ! Diag03182_if1
!
! ----------------------------------------------------------------------
! DIAG.03183 Copy qCL INC: bdy layer + PC2 condensation to stashwork
! ----------------------------------------------------------------------
      item = 183
! Diag03183_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Diag03183_do1:
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              l = i + (j-1)*row_length + (k-1)*rows*row_length
              interp_data_3(l)  = qcl_incr_diagnostic(i, j, k)
            End Do
          End Do
        End Do  ! Diag03183_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       interp_data_3,                                             &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
        cmessage="diag_bl  : error in copydiag_3d(qCL INC: bdy lay+PC2)"
        End if
      End if  ! Diag03183_if1
!
! ----------------------------------------------------------------------
! DIAG.03184 Copy qCF INC: bdy layer to stashwork
! ----------------------------------------------------------------------
      item = 184
! Diag03184_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       qcf_incr_diagnostic,                                       &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="diag_bl  : error in copydiag_3d(Qcf INC: bdy layer)"
        End if
      End if  ! Diag03184_if1
!
! ----------------------------------------------------------------------
! DIAG.03185 Copy U Wind INC: bdy layer to stashwork
! ----------------------------------------------------------------------
      item = 185  ! u wind increment
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        u_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 185)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 186  ! v wind increment
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        v_incr_diagnostic,                                        &
     &        row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 186)"//cmessage
         End if

      Endif  !  sf(item,sect)
!
      item = 192  ! total cloud fraction increment
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cf_incr_diagnostic,                                       &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 192)"//cmessage
         End if

      Endif  !  sf(item,sect)

       item = 193  ! liquid cloud fraction increment
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cfl_incr_diagnostic,                                      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 193)"//cmessage
         End if

      Endif  !  sf(item,sect)
       item = 194  ! ice cloud fraction increment
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cff_incr_diagnostic,                                      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 194)"//cmessage
         End if

      Endif  !  sf(item,sect)
!
! ----------------------------------------------------------------------
! DIAG.03189 Copy Tl(Liq) INC: bdy layer to stashwork
! ----------------------------------------------------------------------
      item = 189
! Diag03189_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Diag03189_do1:
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              l = i + (j-1)*row_length + (k-1)*rows*row_length

              interp_data_3(l) = T_incr_diagnostic(i, j, k)             &
     &                         - LCRCP * qcl_incr_diagnostic(i, j, k)
            End Do
          End Do
        End Do  ! Diag03189_do1
!
! Diag03189_if2:
        If (model_levels   >    wet_model_levels)  Then
!
! Diag03189_do2:
          Do k = (wet_model_levels + 1), model_levels
            Do j = 1, rows
              Do i = 1, row_length
                l = i + (j-1)*row_length + (k-1)*rows*row_length
                interp_data_3(l) = T_incr_diagnostic(i, j, k)
              End Do
            End Do
          End Do  ! Diag03189_do2
!
        End if  ! Diag03189_if2
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       interp_data_3,                                             &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
        cmessage="diag_bl  : error in copydiag_3d(Tl(Liq) INC: bdy lay)"
        End if
      End if  ! Diag03189_if1
!
! ----------------------------------------------------------------------
! DIAG.03190 Copy qT(Liq) INC: bdy layer to stashwork
! ----------------------------------------------------------------------
      item = 190
! Diag03190_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Diag03190_do1:
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              l = i + (j-1)*row_length + (k-1)*rows*row_length
              interp_data_3(l)  =   q_incr_diagnostic(i, j, k)          &
     &                          +   qcl_incr_diagnostic(i, j, k)
            End Do
          End Do
        End Do  ! Diag03190_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       interp_data_3,                                             &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
        cmessage="diag_bl  : error in copydiag_3d(QT(Liq) INC: bdy lay)"
        End if
      End if  ! Diag03190_if1
!
! ----------------------------------------------------------------------
!  10m x Wind
! ----------------------------------------------------------------------
! Item 225 u10m

      If (icode <= 0 .and. sf(225,3)) Then

! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: uc_to_ub
          CALL  uC_to_uB(u10m,                                          &
     &                   row_length,rows,n_rows,1,off_x,off_y,          &
     &                   STASHwork(si(225,3,im_index)))

      End if

      If (sf(209,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(209,3,im_index)),u10m,              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,209,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(u10m) "
         Endif

      End if ! sf(209,3)


! ----------------------------------------------------------------------
!  10m y Wind
! ----------------------------------------------------------------------
! Item 226 v10m

      If (icode <= 0 .and. sf(226,3)) Then

! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: vc_to_vb
          CALL  vC_to_vB(v10m,                                          &
     &                   row_length,n_rows,1,off_x,off_y,               &
     &                   STASHwork(si(226,3,im_index)))

      End if

      If (sf(210,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(210,3,im_index)),v10m,              &
     &        row_length,n_rows,0,0,0,0, at_extremity,                  &
     &        atmos_im,3,210,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(v10m) "
         Endif

      End if ! sf(210,3)


! ----------------------------------------------------------------------
!  10m Wind Speed on B-grid
! ----------------------------------------------------------------------
! Item 227

      If (icode <= 0 .and. sf(227,3)) Then

         ! Horizontal interpolation from 'C' to 'B' grid
! DEPENDS ON: uc_to_ub
         CALL uC_to_uB(u10m,                                            &
     &        row_length,rows,n_rows,1,off_x,off_y,                     &
     &        u10mB)
! DEPENDS ON: vc_to_vb
         CALL vC_to_vB(v10m,                                            &
     &        row_length,n_rows,1,off_x,off_y,                          &
     &        v10mB)

         ! Calculate wind speed
         Do j = 1,n_rows
           Do i = 1,row_length
             ws10mB(i,j) = sqrt(u10mB(i,j)*u10mB(i,j)                   &
     &                         +v10mB(i,j)*v10mB(i,j))
           Enddo !i
         Enddo !j

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(227,3,im_index)), ws10mB,           &
     &        row_length,n_rows,0,0,0,0, at_extremity,                  &
     &        atmos_im,3,227,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(10m wnd spd B) "
         Endif

      End if ! sf(227,3)

! ----------------------------------------------------------------------
!  10m Wind Speed on C-grid P points for ocean
! ----------------------------------------------------------------------
! Item 230

      If (icode <= 0 .and. sf(230,3)) Then

      Allocate(u10m_p(row_length,rows))
      Allocate(v10m_p(row_length,rows))

! Set allocated arrays to zero to make sure any data left in the memory
! allocated to them is purged
      u10m_p(:,:)   =0.0
      v10m_p(:,:)   =0.0

      Allocate(u10m_halo(1-off_x:row_length+off_x, 1-off_y:rows+off_y))
      Allocate(v10m_halo(1-off_x:row_length+off_x,                      &
     &                                           1-off_y:n_rows+off_y))

! Set allocated arrays to zero to make sure any data left in the memory
! allocated to them is purged
      u10m_halo(:,:)=0.0
      v10m_halo(:,:)=0.0

      Do y=1,rows
        Do z=1,row_length
          u10m_halo(z,y)=u10m(z,y)
        End Do
      End Do
      Do y=1,n_rows
        Do z=1,row_length
          v10m_halo(z,y)=v10m(z,y)
        End Do
      End Do

      ! Update halos for u10m_halo
! DEPENDS ON: swap_bounds
        Call Swap_bounds(u10m_halo,row_length,rows,1,                   &
     &                 off_x,off_y,fld_type_u,.true.)
      ! Update halos for v10m_halo
! DEPENDS ON: swap_bounds
        Call Swap_bounds(v10m_halo,row_length,n_rows,1,                 &
     &                 off_x,off_y,fld_type_v,.true.)

! interpolate u and v to p grid.
! DEPENDS ON: u_to_p
      Call U_TO_P(u10m_halo,row_length,rows,1,                          &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, u10m_p)

! DEPENDS ON: v_to_p
      Call V_TO_P(v10m_halo,row_length,rows,n_rows,1,                   &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, v10m_p)

      If (MODEL_DOMAIN  ==  mt_global) Then
! Overwrite values of U_P, V_P at the poles with the magnitude of
! the vector wind.
! DEPENDS ON: polar_vector_wind_n
        Call Polar_vector_wind_n(                                       &
     &                       v10m_halo,                                 &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, 1 , mag_vector_np,                 &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       off_x, off_y, global_row_length,           &
     &                       proc_row_group, at_extremity)

        If (at_extremity(PSouth) ) Then
            Do I=1,ROW_LENGTH
              v10m_p(I,1) = MAG_VECTOR_SP(1)
              u10m_p(I,1) = 0.0
            End Do
        End If

        If (at_extremity(PNorth) ) Then
          Do I=1,ROW_LENGTH
            v10m_p(I,rows) = MAG_VECTOR_NP(1)
            u10m_p(I,rows) = 0.0
          End Do
        End If

      End If

      Deallocate(u10m_halo)
      Deallocate(v10m_halo)

! Calculate 10m wind speed at p C-grid points
       Do j = 1,rows
          Do i = 1,row_length
            ws10m_p(i,j) = sqrt(u10m_p(i,j)*u10m_p(i,j)                 &
     &                         +v10m_p(i,j)*v10m_p(i,j))
          End Do !i
        End Do !j
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(230,3,im_index)),ws10m_p,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,230,                                           &
     &        icode,cmessage)
         If (icode > 0) Then
            cmessage=": error in copydiag(ws10m_p C-grid  ) "
         Endif

         Deallocate(u10m_p)
         Deallocate(v10m_p)

      End If ! sf(230,3)

! ----------------------------------------------------------------------
!   1.5m liquid temperature
! ----------------------------------------------------------------------

! 1.5m liquid temperature
      item =   254  ! tl at 1.5m
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,sect,im_index)),               &
     &        t1p5m,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 254)"//cmessage
         Endif

      End if

! ----------------------------------------------------------------------
!   1.5m  total water
! ----------------------------------------------------------------------

! 1.5m total water
      item =   255  ! qt at 1.5m
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,sect,im_index)),               &
     &        q1p5m,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 255)"//cmessage
         Endif

      End if


! ----------------------------------------------------------------------
!  conversion of 1.5m liquid temperature, total water
! ----------------------------------------------------------------------
      If (icode <= 0 .and. sq_T1p5) Then
!
! Use the diagnostic cloud scheme
! Convert 1.5m TL and QT to T and q, assuming p_star and rhcrit at
!  model level 1 approximates to 1.5m.

! DEPENDS ON: ls_cld
        CALL ls_cld(                                                    &
     &   p_star, rhcpt,  1, boundary_layer_levels, row_length, rows     &
     &  ,rhc_row_length, rhc_rows                                       &
     &  ,cloud_fraction_method,overlap_ice_liquid                       &
     &  ,ice_fraction_method,ctt_weight,t_weight                        &
     &  ,qsat_fixed,sub_cld                                             &
     &  ,ntml, cumulus, l_eacf                                          &
                                                 ! in
     &  ,l_mixing_ratio                                                 &
                                                 ! in
     &  ,t1p5m                                                          &
                                                 ! in/out
     &  ,interp_data_3                                                  &
                                                 ! work_cf
     &  ,q1p5m                                                          &
                                                 ! in/out
     &  ,qcf                                                            &
                                                 ! in
     &  ,qcl1p5m                                                        &
     &  ,interp_data_3(2*row_length*rows+1)                             &
                                                 ! work_cfl
     &  ,interp_data_3(3*row_length*rows+1)                             &
                                                 ! work_cff
     &  ,icode )

      Endif ! on 1.5m STASHflags


! ----------------------------------------------------------------------
!  1.5m Relative humidity
! ----------------------------------------------------------------------
! Item 245: Calculate relative humidity at 1.5m from q1p5m and t1p5m
!   Note that surface pressure is used to determine saturation
!   humidity (instead of p at 1.5m), but this is an extremely close
!   approximation.

      item =   245  ! relative humidity at 1.5m
      If (icode <= 0 .and. sf(item,sect)) Then

!  Find humidity saturation at 1.5m, store in interp_data work.
!  Q1.5m is always specific humidity (conversion in IMP_SOLVER)
!  so we want the specific qsat
! DEPENDS ON: qsat_mix
         Call QSAT_mix(interp_data,t1p5m,p_star,row_length*rows         &
     &                 ,.false.)

         Do j = 1,rows
         Do i = 1,row_length

!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
            interp_data(i,j) = Max( 0. , q1p5m(i,j) ) * 100.            &
     &                                              / interp_data(i,j)

         Enddo ! i
         Enddo ! j

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,sect,im_index)),               &
     &        interp_data,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 245)"//cmessage
         Endif

      Endif  !  sf(item,sect)


! ----------------------------------------------------------------------
!  1.5m Temperature
! ----------------------------------------------------------------------
! Item 236 T1p5m

      If (icode <= 0 .and. sf(236,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(236,3,im_index)),T1p5m,             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,236,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 236)"
         Endif

      End if
! Obukhov length
      item=464
      If (icode <= 0 .and. sf(item,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,3,im_index)),BL_diag%oblen,    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 464)"
         Endif

      End if

! Friction velocity
      item=465
      If (icode <= 0 .and. sf(item,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,3,im_index)),BL_diag%ustar,    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 465)"
         Endif

      End if

! Surface buoyancy flux
      item=467
      If (icode <= 0 .and. sf(item,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,3,im_index)),BL_diag%wbsurf,   &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 467)"
         Endif

      End if

! Gradient Richardson number
      item =   468
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        BL_diag%gradrich,                                         &
     &        row_length,rows,boundary_layer_levels,                    &
     &        0,0,0,0, at_extremity,                                    &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 468)"//cmessage
        End if

      Endif  !  sf(item,sect)

! Convective velocity scale
      item=466
      If (icode <= 0 .and. sf(item,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,3,im_index)),BL_diag%wstar,    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 466)"
         Endif

      End if

! Stratification
      item =   469
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        BL_diag%dbdz,                                             &
     &        row_length,rows,boundary_layer_levels,                    &
     &        0,0,0,0, at_extremity,                                    &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 469)"//cmessage
        End if

      Endif  !  sf(item,sect)

! Modulus of shear
      item =   470
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        BL_diag%dvdzm,                                            &
     &        row_length,rows,boundary_layer_levels,                    &
     &        0,0,0,0, at_extremity,                                    &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 470)"//cmessage
        End if

      Endif  !  sf(item,sect)

! Momentum diffusivity
      item =   471
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        BL_diag%rhokm,                                            &
     &        row_length,rows,boundary_layer_levels,                    &
     &        0,0,0,0, at_extremity,                                    &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 471)"//cmessage
        End if

      Endif  !  sf(item,sect)

! Heat diffusivity
      item =   472
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        BL_diag%rhokh,                                            &
     &        row_length,rows,boundary_layer_levels,                    &
     &        0,0,0,0, at_extremity,                                    &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 472)"//cmessage
        End if

      Endif  !  sf(item,sect)

! Turbulent kinetic energy
      item =   473
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        BL_diag%tke,                                              &
     &        row_length,rows,boundary_layer_levels,                    &
     &        0,0,0,0, at_extremity,                                    &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 473)"//cmessage
        End if

      Endif  !  sf(item,sect)

!     Orographic stress (x component)
      item =   474
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        BL_diag%ostressx,                                         &
     &        row_length,rows,boundary_layer_levels,                    &
     &        0,0,0,0, at_extremity,                                    &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 474)"//cmessage
        End if

      Endif  !  sf(item,sect)

!     Orographic stress (y component)
      item =   475
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        BL_diag%ostressy,                                         &
     &        row_length,rows,boundary_layer_levels,                    &
     &        0,0,0,0, at_extremity,                                    &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 475)"//cmessage
        End if

      Endif  !  sf(item,sect)

! Height of decoupled layer base (=SML top if doesn't exist)
      item=360
      If (icode <= 0 .and. sf(item,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,3,im_index)),BL_diag%dscbase,  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 360)"
         Endif

      End if

! Height of stratocumulus cloud base
      item=361
      If (icode <= 0 .and. sf(item,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,3,im_index)),BL_diag%cldbase,  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 361)"
         Endif

      End if

! Parametrized entrainment rate for surface-based mixed layer
      item=362
      If (icode <= 0 .and. sf(item,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,3,im_index)),BL_diag%weparm,   &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 362)"
         Endif

      End if

! Parametrized entrainment rate for decoupled mixed layer, SML if not
      item=363
      If (icode <= 0 .and. sf(item,3)) Then

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(item,3,im_index)),BL_diag%weparm_dsc,&
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 363)"
         Endif

      End if



! Gust Diagnostic Calculation
      item=463
      if (icode <= 0 .and. sf(item,3)) then

      allocate(gust_wind(row_length,rows))
      allocate(std_dev(row_length,rows))
      allocate(u10m_p(row_length,rows))
      allocate(v10m_p(row_length,rows))

! Set allocated arrays to zero to make sure any data left in the memory
! allocated to them is purged
      gust_wind(:,:)=0.0
      std_dev(:,:)  =0.0
      u10m_p(:,:)   =0.0
      v10m_p(:,:)   =0.0

      do i_g=1,row_length
        do j_g=1,rows
         if (bl_diag%oblen(i_g,j_g) > 0)then
           std_dev(i_g,j_g)=gust_const*bl_diag%ustar(i_g,j_g)
         elseif (bl_diag%oblen(i_g,j_g) < 0)then
           std_dev(i_g,j_g)=gust_const*bl_diag%ustar(i_g,j_g)* (        &
     &      (  1 - ml_depth(i_g,j_g)/(24*bl_diag%oblen(i_g,j_g)) )      &
     &       **(1.0/3.0) )
         elseif(bl_diag%oblen(i_g,j_g) == 0)then
           write(6,*)'error obukhov length is zero'
         endif
        enddo
      enddo

      allocate(u10m_halo(1-off_x:row_length+off_x, 1-off_y:rows+off_y))
      allocate(v10m_halo(1-off_x:row_length+off_x,                      &
     &                                           1-off_y:n_rows+off_y))

! Set allocated arrays to zero to make sure any data left in the memory
! allocated to them is purged
      u10m_halo(:,:)=0.0
      v10m_halo(:,:)=0.0

      do y=1,rows
        do z=1,row_length
          u10m_halo(z,y)=u10m(z,y)
        enddo
      enddo
      do y=1,n_rows
        do z=1,row_length
          v10m_halo(z,y)=v10m(z,y)
        enddo
      enddo

      ! Update halos for u10m_halo
! DEPENDS ON: swap_bounds
        CALL Swap_bounds(u10m_halo,row_length,rows,1,                   &
     &                 off_x,off_y,fld_type_u,.true.)
      ! Update halos for v10m_halo
! DEPENDS ON: swap_bounds
        CALL Swap_bounds(v10m_halo,row_length,n_rows,1,                 &
     &                 off_x,off_y,fld_type_v,.true.)

! interpolate u and v to p grid.
! DEPENDS ON: u_to_p
      Call U_TO_P(u10m_halo,row_length,rows,1,                          &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, u10m_p)


! DEPENDS ON: v_to_p
      Call V_TO_P(v10m_halo,row_length,rows,n_rows,1,                   &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, v10m_p)

      IF(MODEL_DOMAIN  ==  mt_global) THEN
! Overwrite values of U_P, V_P at the poles with the magnitude of
! the vector wind.
! DEPENDS ON: polar_vector_wind_n
        Call Polar_vector_wind_n(                                       &
     &                       v10m_halo,                                 &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, 1 , mag_vector_np,                 &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       off_x, off_y, global_row_length,           &
     &                       proc_row_group, at_extremity)

        If (at_extremity(PSouth) ) Then
            DO I=1,ROW_LENGTH
              v10m_p(I,1) = MAG_VECTOR_SP(1)
              u10m_p(I,1) = 0.0
            END DO
        End If

        If (at_extremity(PNorth) ) Then
            DO I=1,ROW_LENGTH
              v10m_p(I,rows) = MAG_VECTOR_NP(1)
              u10m_p(I,rows) = 0.0
            END DO
        End If

      ENDIF




      deallocate(u10m_halo)
      deallocate(v10m_halo)

        Do j_g = 1,rows
          Do i_g = 1,row_length
!Calculate the scaling factor
            scale=(1/vkman)*log( (5*exp(vkman*c_ugn)+z0m_eff(i_g,j_g))/ &
     &                        (5+z0m_eff(i_g,j_g)))

            gust_wind(i_g,j_g) = (sqrt(u10m_p(i_g,j_g)*u10m_p(i_g,j_g)  &
     &                         + v10m_p(i_g,j_g)* v10m_p(i_g,j_g))) +   &
     &                           scale*std_dev(i_g,j_g)
          Enddo !i_g
        Enddo !j_g

!      if (icode <= 0 .and. sf(item,3)) then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,3,im_index)),gust_wind,        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,item,                                          &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage=": error in copydiag(item 463)"
         Endif

         deallocate(gust_wind)
         deallocate(std_dev)
         deallocate(u10m_p)
         deallocate(v10m_p)

      End if


! Silhouette Orographic Roughness (A/S)
! ----------------------------------------------------------------------
! Item 251

      IF (SF(251,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(251,3,im_index)),                              &
     &       sil_orog_land,                                             &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF


! Orographic Roughness Peak to Trough Height
! ----------------------------------------------------------------------
! Item 252

      IF (SF(252,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(252,3,im_index)),                              &
     &       ho2r2_orog,                                                &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF



! ----------------------------------------------------------------------
!  1.5m Specific Humidity
! ----------------------------------------------------------------------
! Item 237 q1p5m

      If (icode <= 0 .and. sf(237,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(237,3,im_index)),q1p5m,             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,237,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 237)"
         Endif

      End if


! ----------------------------------------------------------------------
!  1.5m Fog Fraction and Mist Fraction
! ----------------------------------------------------------------------
!  Items 248 and 253

      If (icode <= 0 .and. (sf(248,3) .or. sf(253,3)) ) Then
        Do i = 1,row_length
          Do j = 1,rows
            Do k=1,n_vis_thresh
              Vis_Threshold(i,j,1,k)=vis_thresh(k)
            Enddo
          Enddo
        Enddo
! DEPENDS ON: fog_fr
        Call fog_fr (p_star, RHcrit,1,                                  &
     &       row_length * rows,                                         &
     &       T1p5m, Aerosol, L_murk,                                    &
     &       q1p5m, Qcl1p5m, Qcf,                                       &
     &       vis_threshold,PVis,n_vis_thresh)

        If (icode <= 0 .and. sf(248,3)) Then
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(248,3,im_index)), PVIS(1,1,1),     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,248,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 248)"
         Endif
        Endif     ! sf(248,3)

        If (icode <= 0 .and. sf(253,3)) Then
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(253,3,im_index)), PVIS(1,1,2),     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,253,                                           &
     &        icode,cmessage)

          If (icode  >   0) Then
            cmessage=": error in copydiag(item 253)"
          Endif
        Endif    ! sf(253,3)

      End if


! ----------------------------------------------------------------------
!  1.5m Dewpoint
! ----------------------------------------------------------------------
! Item 250

      If (icode <= 0 .and. sf(250,3)) Then

! DEPENDS ON: dewpnt
         Call DEWPNT (q1p5m, p_star, T1p5m,                             &
     &               row_length * rows, interp_data)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(250,3,im_index)), interp_data,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,250,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 250)"
         Endif

      End if


! ----------------------------------------------------------------------
!  Visibility in precipitation diagnostics
! ----------------------------------------------------------------------
! Items 281, 282, 283, 284 and 285.

       If (icode <= 0 .and. (sf(281,3) .or. sf(282,3) .or. sf(283,3)    &
     &                       .or. sf(284,3) .or. sf(285,3)) ) Then
!
!     Calculate scattering coefficients due to precipitation
!
! DEPENDS ON: beta_precip
        Call Beta_precip(ls_Rain, ls_Snow                               &
     &                 ,Conv_Rain, Conv_Snow, qcf(1,1,1)                &
     &                 ,rho1, T1p5m, p_star                             &
     &                 ,plsp,cca_2d,.false.,.true.                      &
     &                 ,row_length*rows,row_length*rows,1               &
     &                 ,x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic        &
     &                 ,lsp_ei,lsp_fi,lsp_eic,lsp_fic                   &
     &                 ,Beta_ls_Rain, Beta_ls_Snow                      &
     &                 ,Beta_C_Rain, Beta_C_Snow, icode)

      Endif

      If (icode <= 0 .and. (sf(282,3) .or. sf(283,3)) ) Then
!
!     Calculate screen level probability of vis less than thresholds
!
! DEPENDS ON: calc_vis_prob
        Call Calc_vis_prob(p_star,                                      &
     &       rhcrit,1,                                                  &
     &       row_length*rows,row_length*rows,                           &
     &       T1p5m,aerosol,l_murk,                                      &
     &       q1p5m,qcl1p5m,qcf,                                         &
     &       vis_thresh,n_vis_thresh,                                   &
     &       plsp,cca_2d,.false.,                                       &
     &       beta_ls_rain,beta_ls_snow,                                 &
     &       beta_c_rain,beta_c_snow,                                   &
     &       vis_pr,                                                    &
     &       icode)
      Endif
!
! Item 282 Prob vis<1000 m
!
      If(icode <= 0 .and. sf(282,3)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(282,3,im_index)), vis_pr(1,1,1,1),  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,282,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 282)"
         Endif
      Endif
!
! Item 283 Prob vis<5000 m
!
      If(icode <= 0 .and. sf(283,3)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(283,3,im_index)), vis_pr(1,1,1,2),  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,283,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 283)"
         Endif
      Endif

      If (icode <= 0 .and.                                              &
     &    (sf(281,3) .OR. sf(284,3) .OR. sf(285,3)) ) Then
!
!     Visibility at 1.5 m including precipitation
!
! DEPENDS ON: visbty
         Call VISBTY(                                                   &
     &     p_star, T1p5m, q1p5m, Qcl1p5m, Qcf                           &
                                                      !INPUT
     &     ,Aerosol, calc_prob_of_vis, RHcrit, L_murk                   &
                                                      !INPUT
     &     ,row_length * rows                                           &
                                                      !INPUT
     &     ,Vis_no_precip)                            !OUTPUT
! DEPENDS ON: vis_precip
         Call VIS_PRECIP(Vis_No_Precip                                  &
     &               ,plsp,cca_2d,.false.                               &
     &               ,Beta_ls_Rain, Beta_ls_Snow                        &
     &               ,Beta_C_Rain, Beta_C_Snow                          &
     &               ,row_length*rows,row_length*rows,1                 &
     &               ,vis,Vis_ls_Precip,Vis_C_Precip                    &
     &               ,icode)
      Endif

!
! Item 281 Median Visibility including precip
!
      If(icode <= 0 .and. sf(281,3)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(281,3,im_index)), vis,              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,281,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 281)"
         Endif
      Endif
!
! Item 284 Visibility in ls precip
!
      If(icode <= 0 .and. sf(284,3)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(284,3,im_index)), vis_ls_Precip,    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,284,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 284)"
         Endif
      Endif
!
! Item 285 Visibility in Conv precip
!
      If(icode <= 0 .and. sf(285,3)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(285,3,im_index)), vis_C_Precip,     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,285,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 285)"
         Endif

      End if

! ----------------------------------------------------------------------
!   Sensible Heat flux: surface and model levels
! ----------------------------------------------------------------------
! Items 216,217 array ftl combines surface and model level values
!         (,,1) is the surface value.
!         Since model level fields are on rho_levels and the surface
!         is a theta_level, these are output as 2 separate diagnostics.

      item =   216  ! T_L flux profile on rho_levels
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        ftl,                                                      &
     &        row_length,rows,boundary_layer_levels,                    &
     &        0,0,0,0, at_extremity,                                    &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 216)"//cmessage
        End if

      Endif  !  sf(item,sect)

! surface = ftl(,,1)
      item =   217  ! T_L flux at surface
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,sect,im_index)),               &
     &        ftl,                                                      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 217)"//cmessage
         Endif

      End if

! ----------------------------------------------------------------------
!   Sensible Moisture flux: surface and model levels
! ----------------------------------------------------------------------
! Items 222,223,241 array fqt combines surface and model level values
!         (,,1) is the surface value.
!         Since model level fields are on rho_levels and the surface
!         is a theta_level, these are output as separate diagnostics.

      item =   222  ! q_T flux profile on rho_levels
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        fqt,                                                      &
     &        row_length,rows,boundary_layer_levels,                    &
     &        0,0,0,0, at_extremity,                                    &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 222)"//cmessage
        End if

      Endif  !  sf(item,sect)

! surface = fqt(,,1)
      item =   223  ! T_L flux at surface
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,sect,im_index)),               &
     &        fqt,                                                      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 223)"//cmessage
         Endif

      Endif  !  sf(item,sect)

      item =   241  ! total surface moisture flux per timestep
      If (icode <= 0 .and. sf(item,sect)) Then

! Convert from rate to timestep accumulation explicitly
         Do j = 1,rows
           Do i = 1,row_length
             interp_data(i,j) = fqt(i,j,1) * timestep
           End Do
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(item,sect,im_index)),               &
     &        interp_data,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 241)"//cmessage
         Endif

      Endif  !  sf(item,sect)

! ----------------------------------------------------------------------
!   Latent Heat flux.
! ----------------------------------------------------------------------
! Item 234 surf_latent_heat_flux

      If (icode <= 0 .and. sf(234,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(234,3,im_index)),                   &
     &        surf_latent_heat_flux,                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,234,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="error in copydiag(item 234)"
         Endif

      End if


! ----------------------------------------------------------------------
!  U component of wind stress.
! ----------------------------------------------------------------------
! Item 219 taux

      If (icode <= 0 .and. sf(219,3)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(219,3,im_index)),taux,           &
     &        row_length,rows,boundary_layer_levels,0,0,0,0,            &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,219,3,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,3,219,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage=": error in copydiag_3d(item 219)"
         Endif

      End if


! ----------------------------------------------------------------------
!  V component of wind stress.
! ----------------------------------------------------------------------
! Item 220 tauy

      If (icode <= 0 .and. sf(220,3)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(220,3,im_index)),tauy,           &
     &        row_length,n_rows,boundary_layer_levels,0,0,0,0,          &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,220,3,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,3,220,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage=": error in copydiag_3d(item 220)"
         Endif

      End if


! ----------------------------------------------------------------------
!  Magnitude of wind stress on B-grid
! ----------------------------------------------------------------------
! Item 221 tauB

      If (icode <= 0 .and. sf(221,3)) Then

         ! Horizontal interpolation from 'C' to 'B' grid
! DEPENDS ON: uc_to_ub
         CALL uC_to_uB(taux,                                            &
     &        row_length,rows,n_rows,boundary_layer_levels,             &
     &        off_x,off_y,tauxB)
! DEPENDS ON: vc_to_vb
         CALL vC_to_vB(tauy,                                            &
     &        row_length,n_rows,boundary_layer_levels,                  &
     &        off_x,off_y,tauyB)

         ! Calculate stress magnitude
         Do k = 1,boundary_layer_levels
           Do j = 1,n_rows
             Do i = 1,row_length
               tauB(i,j,k) = sqrt(tauxB(i,j,k)*tauxB(i,j,k)             &
     &                           +tauyB(i,j,k)*tauyB(i,j,k))
             Enddo !i
           Enddo !j
         Enddo !k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(221,3,im_index)),tauB,           &
     &        row_length,n_rows,boundary_layer_levels,0,0,0,0,          &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,221,3,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,3,221,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag_3d(item 221)"
         Endif

      End if

! ----------------------------------------------------------------------
! Section 2.2.9  Soil temperature
! ----------------------------------------------------------------------
! Item 238  t_soil

      If (icode <= 0 .and. sf(238,3)) Then

         Do k = 1, dsm_levels
            Do j= 1, rows
               Do i = 1, row_length
                 l = i + (j-1)*row_length + (k-1)*rows*row_length
                 interp_data_3(l) = rmdi
               End Do
            End Do
            Do i = 1, land_points
               l = land_index(i) + (k-1)*rows*row_length
               interp_data_3(l) = t_soil(i,k)
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(238,3,im_index)),interp_data_3,  &
     &        row_length,rows,dsm_levels,0,0,0,0,                       &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,238,3,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,3,238,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 238)"
         Endif

      End if


! ----------------------------------------------------------------------
!  Mixed Layer depth.
! ----------------------------------------------------------------------
! Item 025 ml_depth

      If (icode <= 0 .and. sf(025,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(025,3,im_index)),                   &
     &        ml_depth,                                                 &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,025,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 025)"
         Endif

      End if


! ----------------------------------------------------------------------
!   wind mixing energy.
! ----------------------------------------------------------------------
! Item 224 wind_mixing_energy

      If (icode <= 0 .and. sf(224,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(224,3,im_index)),                   &
     &        wind_mixing_energy,                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,224,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 224)"
         Endif

      End if


! ----------------------------------------------------------------------
!  Boundary Layer top.
! ----------------------------------------------------------------------
! Item 304 bl_top

      If (icode <= 0 .and. sf(304,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(304,3,im_index)),                   &
     &        bl_top,                                                   &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,304,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 304)"
         Endif

      End if


! ----------------------------------------------------------------------
!  Boundary Layer Type (all 6)
! ----------------------------------------------------------------------

! Item 305...310 bl_type  1..6

      If (icode <= 0 .and. sf(PP_code_bl_types(1),3)) Then

! DEPENDS ON: copydiag
         Call copydiag(                                                 &
     &        STASHwork(si(PP_code_bl_types(1),3,im_index)),            &
     &        bl_type_1,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,PP_code_bl_types(1),                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(bl_type_1)"
         Endif

      End if

      If (icode <= 0 .and. sf(PP_code_bl_types(2),3)) Then

! DEPENDS ON: copydiag
         Call copydiag(                                                 &
     &        STASHwork(si(PP_code_bl_types(2),3,im_index)),            &
     &        bl_type_2,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,PP_code_bl_types(2),                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(bl_type_2)"
         Endif

      End if

      If (icode <= 0 .and. sf(PP_code_bl_types(3),3)) Then

! DEPENDS ON: copydiag
         Call copydiag(                                                 &
     &        STASHwork(si(PP_code_bl_types(3),3,im_index)),            &
     &        bl_type_3,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,PP_code_bl_types(3),                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(bl_type_3)"
         Endif

      End if

      If (icode <= 0 .and. sf(PP_code_bl_types(4),3)) Then

! DEPENDS ON: copydiag
         Call copydiag(                                                 &
     &        STASHwork(si(PP_code_bl_types(4),3,im_index)),            &
     &        bl_type_4,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,PP_code_bl_types(4),                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(bl_type_4)"
         Endif

      End if

      If (icode <= 0 .and. sf(PP_code_bl_types(5),3)) Then

! DEPENDS ON: copydiag
         Call copydiag(                                                 &
     &        STASHwork(si(PP_code_bl_types(5),3,im_index)),            &
     &        bl_type_5,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,PP_code_bl_types(5),                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(bl_type_5)"
         Endif

      End if

      If (icode <= 0 .and. sf(PP_code_bl_types(6),3)) Then

! DEPENDS ON: copydiag
         Call copydiag(                                                 &
     &        STASHwork(si(PP_code_bl_types(6),3,im_index)),            &
     &        bl_type_6,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,PP_code_bl_types(6),                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(bl_type_6)"
         Endif

      End if

      If (icode <= 0 .and. sf(PP_code_bl_types(7),3)) Then

! DEPENDS ON: copydiag
         Call copydiag(                                                 &
     &        STASHwork(si(PP_code_bl_types(7),3,im_index)),            &
     &        bl_type_7,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,PP_code_bl_types(7),                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(bl_type_7)"
         Endif

      End if

! Assign combined boundary layer type diagnostic

      item=476
      If (icode <= 0 .and. sf(item,3)) Then

        Do J=1,ROWS
          Do I=1,ROW_LENGTH

          bl_type_comb(i,j)=                                            &
     &    1.0*bl_type_1(i,j) + 2.0*bl_type_2(i,j) + 3.0*bl_type_3(i,j)  &
     &  + 4.0*bl_type_4(i,j) + 5.0*bl_type_5(i,j) + 6.0*bl_type_6(i,j)  &
     &  + 7.0*bl_type_7(i,j)

          EndDo
        EndDo

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(item,3,im_index)),bl_type_comb,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,item,                                          &
     &        icode,cmessage)

        If (icode  >   0) Then
           cmessage=": error in copydiag(item 476)"
        Endif

      End if  ! item 476

! ----------------------------------------------------------------------
!   Effective roughness length for momentum
! ----------------------------------------------------------------------
! Item 026 z0m_eff

      If (icode <= 0 .and. sf(026,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(026,3,im_index)),                   &
     &        z0m_eff,                                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,026,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 026)"
         Endif

      End if

! ----------------------------------------------------------------------
!   Effective roughness length for heat
! ----------------------------------------------------------------------
! Item 027 z0h_eff

      If (icode <= 0 .and. sf(027,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(027,3,im_index)),                   &
     &        z0h_eff,                                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,027,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 027)"
         Endif

      End if

! ----------------------------------------------------------------------
!   Grid-box mean Vegetative roughness length for momentum
! ----------------------------------------------------------------------
! Item 028 z0m_gb

      If (icode <= 0 .and. sf(028,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(028,3,im_index)),                   &
     &        z0m_gb,                                                   &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,028,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 028)"
         Endif

      End if

! ----------------------------------------------------------------------
!   RIB - measure of surface stability
! ----------------------------------------------------------------------
! Item 208 rib

      If (icode <= 0 .and. sf(208,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(208,3,im_index)),                   &
     &        rib,                                                      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,208,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 208)"
         Endif

      End if

! ----------------------------------------------------------------------
! visibility diagnostics
! ----------------------------------------------------------------------
! Item 247 visibility

      If (icode <= 0 .and. sf(247,3)) Then

! DEPENDS ON: visbty
         Call VISBTY(                                                   &
     &     p_star, T1p5m, q1p5m, Qcl1p5m, Qcf                           &
                                                     !INPUT
     &     ,Aerosol, calc_prob_of_vis, RHcrit, L_murk                   &
                                                     !INPUT
     &     ,row_length * rows                                           &
                                                 !INPUT
     &     ,interp_data)                         !OUTPUT

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(247,3,im_index)),                   &
     &        interp_data,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,247,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 247)"

         Endif

      End if


! Surface evaporation weighted by leads
! ----------------------------------------------------------------------
! Item 232 Surface evaporation weighted by leads

      If (icode <= 0 .and. sf(232,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(232,3,im_index)),                   &
     &        e_sea,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,232,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 232)"
         Endif

      End if


! Sensible heat flux over open sea
! ----------------------------------------------------------------------
! Item 228 Sensible heat flux over open sea

      If (icode <= 0 .and. sf(228,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(228,3,im_index)),                   &
     &        h_sea,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,228,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 228)"
         Endif

      End if

! Evaporation amount from soil surface
! ----------------------------------------------------------------------
! Item 229 soil evaporation

      If (icode <= 0 .and. sf(229,3)) Then

! Convert from rate to timestep accumulation explicitly
!
         Do j = 1,rows
           Do i = 1,row_length
             interp_data(i,j) = es(i,j) * timestep
           End Do
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(229,3,im_index)),                   &
     &        interp_data,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,229,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 229)"
         Endif

      End if

! Sublimation ( if sea-ice included in ocean)
! ----------------------------------------------------------------------
! Item 231 Sublimation

      If (icode <= 0 .and. sf(231,3)) Then

! Convert from rate to timestep accumulation explicitly
         Do j = 1,rows
           Do i = 1,row_length
             interp_data(i,j) = ei(i,j) * timestep
           End Do
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(231,3,im_index)),                   &
     &        interp_data,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,231,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 231)"
         Endif

      End if

! Heat flux through sea-ice ( if sea-ice included in ocean)
! ----------------------------------------------------------------------
! Item 201 Melting

      If (icode <= 0 .and. sf(201,3)) Then
         SEA_ICE_HTF_GB(:,:)=0.0
         DO J=1,ROWS
           DO I=1,ROW_LENGTH
             DO K=1,NICE
               SEA_ICE_HTF_GB(I,J)=SEA_ICE_HTF_GB(I,J)+                 &
     &                                     SEA_ICE_HTF(I,J,K)
             ENDDO
           ENDDO
         ENDDO

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(201,3,im_index)),                   &
     &        SEA_ICE_HTF_GB,                                           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,201,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 201)"
         Endif

      End if


! ----------------------------------------------------------------------
! 3d ice catagory heat flux through sea ice
! ----------------------------------------------------------------------
! Item 256  sea_ice_htf

      If (icode <= 0 .and. sf(256,3)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NICE,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,256,3,im_index)),                       &
     &       PLLNICE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 256 = npp_ft)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NICE
          IF (PLLNICE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: copydiag
            Call copydiag(STASHwork(si(256,3,im_index)+(PSLEVEL_OUT-1)  &
     &        *row_length*rows),                                        &
     &        SEA_ICE_HTF(1,1,PSLEVEL),                                 &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,256,                                           &
     &        icode,cmessage)
            If (icode  >   0) Then
                cmessage=": error in copydiag(item 256)"
            Endif
          ENDIF
        ENDDO
      Endif


! Melting of top of sea-ice ( if sea-ice included in ocean)
! ----------------------------------------------------------------------
! Item 235 Melting

      If (icode <= 0 .and. sf(235,3)) Then
         SICE_MLT_HTF_GB(:,:)=0.0
         DO J=1,ROWS
           DO I=1,ROW_LENGTH
             DO K=1,NICE
               SICE_MLT_HTF_GB(I,J)=SICE_MLT_HTF_GB(I,J)+               &
     &                                  SICE_MLT_HTF(I,J,K)
             ENDDO
           ENDDO
         ENDDO

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(235,3,im_index)),                   &
     &        SICE_MLT_HTF_GB,                                          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,235,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 235)"
         Endif

      End if

! ----------------------------------------------------------------------
! 3d ice catagory sea ice surface melt heat flux
! ----------------------------------------------------------------------
! Item 257  sice_mlt_htf

      If (icode <= 0 .and. sf(257,3)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NICE,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,257,3,im_index)),                       &
     &       PLLNICE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 257 = npp_ft)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NICE
          IF (PLLNICE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: copydiag
            Call copydiag(STASHwork(si(257,3,im_index)+(PSLEVEL_OUT-1)  &
     &        *row_length*rows),                                        &
     &        SICE_MLT_HTF(1,1,PSLEVEL),                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,257,                                           &
     &        icode,cmessage)
            If (icode  >   0) Then
                cmessage=": error in copydiag(item 257)"
            Endif
          ENDIF
        ENDDO
      End if


! Heat flux due to melting of snow (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 258 Melting

      If (icode <= 0 .and. sf(258,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(258,3,im_index)),                   &
     &        snomlt_surf_htf,                                          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,258,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 258)"
         Endif

      End if


! Heat flux from surface to deep soil level 1 (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 202 Surf soil flux

      If (sf(202,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(202,3,im_index)),                   &
     &        surf_ht_flux,                                             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,202,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=                                                   &
     &       "imp_ctl  : error in copydiag(item 202 = surf_ht_flux)"
            goto 9999
         Endif

      End if



! Evaporation from soil surface : rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 296 Soil evap

      If (sf(296,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(296,3,im_index)),                   &
     &        es,                                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,296,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=                                                   &
     &       "imp_ctl  : error in copydiag(item 296 = es)"
            goto 9999
         Endif

      End if



! Evaporation from canopy : rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 297 Canopy evap

      If (sf(297,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(297,3,im_index)),                   &
     &        ecan,                                                     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,297,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=                                                   &
     &       "imp_ctl  : error in copydiag(item 297 = ecan)"
            goto 9999
         Endif

      End if



! Surface sublimation : rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 298 Surf sublim

      If (sf(298,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(298,3,im_index)),                   &
     &        ei,                                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,298,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=                                                   &
     &       "imp_ctl  : error in copydiag(item 298 = ei)"
            goto 9999
         Endif

      End if



! TOA outgoing longwave radiation (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 332 TOA outgoing LW

      If (sf(332,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(332,3,im_index)),                   &
     &        olr,                                                      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,332,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=                                                   &
     &       "imp_ctl  : error in copydiag(item 332 = olr)"
            goto 9999
         Endif

      End if

! Land Mean Potential Evaporation: rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 334 Land potential evap

      If (sf(334,3)) Then

         do l=1,land_points
           epot_land(l)=0.0
           do n=1,ntiles
             epot_land(l)=epot_land(l)+tile_frac(l,n)*epot_tile(l,n)
           enddo
         enddo

! DEPENDS ON: from_land_points
         CALL FROM_LAND_POINTS (                                        &
     &      STASHWORK(SI(334,3,im_index)),                              &
     &       epot_land,                                                 &
     &       land_sea_mask,row_length*rows,land_points_dum)

      End if

! Potential evaporation on tiles
! ----------------------------------------------------------------------
! Item 335 Tiled potential evap

      IF (SF(335,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,335,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 335 = epot_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(335,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),epot_tile(1,PSLEVEL_OUT),            &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF


! Canopy conductance (m/s)
! ----------------------------------------------------------------------
! Item 259 Canopy conductance

      IF (SF(259,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(259,3,im_index)),                              &
     &       gs,                                                        &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF



! Gross primary productivity (Kg C/m2/s)
! ----------------------------------------------------------------------
! Item 261 GPP

      IF (SF(261,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(261,3,im_index)),                              &
     &       gpp,                                                       &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF



! Net primary productivity (Kg C/m2/s)
! ----------------------------------------------------------------------
! Item 262 NPP

      IF (SF(262,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(262,3,im_index)),                              &
     &       npp,                                                       &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF



! Plant respiration (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 263 Plant respiration

      IF (SF(263,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(263,3,im_index)),                              &
     &       resp_p,                                                    &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF



! Soil respiration (Kg C/m2/s)
! ----------------------------------------------------------------------
! Item 293 Soil respiration
! Lestevens March 2010: added cable switch to if loop.
      IF (SF(293,3)) THEN
!        IF (DIM_CS1  ==  4) THEN
        IF (l_cable .or. DIM_CS1  ==  4) THEN
! DEPENDS ON: from_land_points
          CALL FROM_LAND_POINTS (                                       &
     &        STASHWORK(SI(293,3,im_index)),                            &
     &         resp_s_tot,                                              &
     &         land_sea_mask,row_length*rows,land_points_dum)
        ELSE
! DEPENDS ON: from_land_points
          CALL FROM_LAND_POINTS (                                       &
     &        STASHWORK(SI(293,3,im_index)),                            &
     &         resp_s,                                                  &
     &         land_sea_mask,row_length*rows,land_points_dum)
        ENDIF
      END IF

! Soil respiration on tiles (Kg C/m2/s)
! -----------------------------------------------
! Item 173 - kdcorbin, 10/10
     IF (SF(173,3)) THEN
        ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,173,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item 173)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(173,3,im_index)+(PSLEVEL_OUT-1)*row_length*rows) &
      &          , resp_s_tile(1,PSLEVEL_OUT)           &
      &          ,land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
     ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Soil carbon (Kg C m-2)
!--------------------------------------
! ITEM 320: TOTAL SOIL CARBON CONTENT

      IF (SF(320,3)) THEN
        IF (DIM_CS1  ==  4) THEN
          DO i=1,LAND_POINTS
            CS_TOT(I) = CS(i,1) + CS(i,2) + CS(i,3) + CS(i,4)
          ENDDO
! DEPENDS ON: from_land_points
          CALL FROM_LAND_POINTS (                                       &
     &        STASHWORK(SI(320,3,im_index)),                            &
     &         CS_TOT,                                                  &
     &         land_sea_mask,row_length*rows,land_points_dum)
        ELSE
! DEPENDS ON: from_land_points
          CALL FROM_LAND_POINTS (                                       &
     &        STASHWORK(SI(320,3,im_index)),                            &
     &         CS,                                                      &
     &         land_sea_mask,row_length*rows,land_points_dum)
        ENDIF
      END IF

! Soil carbon (Kg C m-2)
!----------------------------------------------------
! ITEMS 477-480: INDIVIDUAL POOL SOIL CARBON CONTENT

! DPM
      IF (SF(477,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(477,3,im_index)),                              &
     &       CS(:,1),                                                   &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF

! RPM
      IF (SF(478,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(478,3,im_index)),                              &
     &       CS(:,2),                                                   &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF

! BIO
      IF (SF(479,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(479,3,im_index)),                              &
     &       CS(:,3),                                                   &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF

!HUM
      IF (SF(480,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(480,3,im_index)),                              &
     &       CS(:,4),                                                   &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF

! Soil respiration (Kg C m-2 s-1)
!-------------------------------------------------------------
! ITEMS 481-484: INDIVIDUAL POOL SOIL RESPIRATION

! DPM
      IF (SF(481,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(481,3,im_index)),                              &
     &       RESP_S(:,1),                                               &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF

! RPM
      IF (SF(482,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(482,3,im_index)),                              &
     &       RESP_S(:,2),                                               &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF

! BIO
      IF (SF(483,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(483,3,im_index)),                              &
     &       RESP_S(:,3),                                               &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF

!HUM
      IF (SF(484,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(484,3,im_index)),                              &
     &       RESP_S(:,4),                                               &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF

! Canopy evaporation on tiles
! ----------------------------------------------------------------------
! Item 287 Tiled canopy evap

      IF (SF(287,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,287,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 287 = ecan_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(287,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),ecan_tile(1,PSLEVEL_OUT),            &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Transpiration + soil evaporation on tiles
! ----------------------------------------------------------------------
! Item 288 Tiled soil evap

      IF (SF(288,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,288,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 288 = esoil_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(288,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),esoil_tile(1,PSLEVEL_OUT),           &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Surface sensible heat flux on tiles (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 290 Tiled surface heat flux

      IF (SF(290,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,290,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 290 = ftl_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(290,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),ftl_tile(1,PSLEVEL_OUT),             &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Bulk Richardson number on tiles
! ----------------------------------------------------------------------
! Item 294 Tiled RIB

      IF (SF(294,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,294,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 294 = rib_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(294,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),rib_tile(1,PSLEVEL_OUT),             &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Surface net radiation on tiles (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 314 Tiled net radiation

      IF (SF(314,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,314,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 314 = radnet_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(314,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),radnet_tile(1,PSLEVEL_OUT),          &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Surface temperature on tiles
! ----------------------------------------------------------------------
! Item 316 Tiled TSTAR

      IF (SF(316,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,316,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 316 = tstar_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(316,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),tstar_tile(1,PSLEVEL_OUT),           &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF


! Surface tile fractions
! ----------------------------------------------------------------------
! Item 317 Surface fractions

      IF (sf(317,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(ntiles,len_stlist,                         &
     &    stlist(1,stindex(1,317,3,im_index)),                          &
     &    plltile,stash_pseudo_levels,num_stash_pseudo,                 &
     &    icode,cmessage)
        IF (icode >  0) THEN
          cmessage=                                                     &
     &    "imp_ctl  : error in set_pseudo_list (item 317 = tile_frac)"
          goto 9999
        END IF
        pslevel_out=0
        DO pslevel=1,ntiles
          IF (plltile(pslevel)) THEN
            pslevel_out = pslevel_out+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS(                                      &
     &        stashwork(si(317,3,im_index)+(pslevel_out-1)              &
     &        *row_length*rows),tile_frac(1,pslevel_out),               &
     &        land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF


! LAI on tiles
! ----------------------------------------------------------------------
! Item 318 LAI on tiles

      IF (sf(318,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(npft,len_stlist,                           &
     &    stlist(1,stindex(1,318,3,im_index)),                          &
     &    pllpft,stash_pseudo_levels,num_stash_pseudo,                  &
     &    icode,cmessage)
        IF (icode >  0) THEN
          cmessage=                                                     &
     &    "imp_ctl  : error in set_pseudo_list (item 318 = lai_ft)"
          goto 9999
        END IF
        pslevel_out=0
        DO pslevel=1,npft
          IF (pllpft(pslevel)) THEN
            pslevel_out = pslevel_out+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS(                                      &
     &        stashwork(si(318,3,im_index)+(pslevel_out-1)              &
     &        *row_length*rows),lai_ft(1,pslevel_out),                  &
     &        land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF


! Canopy height on tiles
! ----------------------------------------------------------------------
! Item 319 Canopy height on tiles

      IF (sf(319,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(npft,len_stlist,                           &
     &    stlist(1,stindex(1,319,3,im_index)),                          &
     &    pllpft,stash_pseudo_levels,num_stash_pseudo,                  &
     &    icode,cmessage)
        IF (icode >  0) THEN
          cmessage=                                                     &
     &    "imp_ctl  : error in set_pseudo_list (item 319 = canht_ft)"
          goto 9999
        END IF
        pslevel_out=0
        DO pslevel=1,npft
          IF (pllpft(pslevel)) THEN
            pslevel_out = pslevel_out+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS(                                      &
     &        stashwork(si(319,3,im_index)+(pslevel_out-1)              &
     &        *row_length*rows),canht_ft(1,pslevel_out),                &
     &        land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF


! "Stomatal" conductance (m/s) on PFTs [pfv: CMIP5 access-1.3 here]
! ----------------------------------------------------------------------
! Item 462 Stomatal conductance on PFTs
! Just consider PFTs as only the stomatal conductance over
! vegetated tiles is needed.
!     IF (sf(462,3)) THEN
! DEPENDS ON: set_pseudo_list
!       CALL SET_PSEUDO_LIST(npft,len_stlist,                           &
!    &    stlist(1,stindex(1,462,3,im_index)),                          &
!    &    pllpft,stash_pseudo_levels,num_stash_pseudo,                  &
!    &    icode,cmessage)
!       IF (icode >  0) THEN
!         cmessage=                                                     &
!    &    "imp_ctl  : error in set_pseudo_list (item 462 = gc)"
!         goto 9999
!       END IF
!      pslevel_out=0
!       DO pslevel=1,npft
!         IF (pllpft(pslevel)) THEN
!           pslevel_out = pslevel_out+1
! DEPENDS ON: from_land_points
!           CALL FROM_LAND_POINTS(                                      &
!    &        stashwork(si(462,3,im_index)+(pslevel_out-1)              &
!    &        *row_length*rows),gc(1,pslevel_out),                      &
!    &        land_sea_mask,row_length*rows,land_points_dum)
!         END IF
!       END DO
!     END IF


! "Stomatal" conductance (m/s) on TILES (pfv: save on all tiles, 21oct13)
! ----------------------------------------------------------------------
! Item 462 Stomatal conductance on TILES
! Just consider PFTs as only the stomatal conductance over
! vegetated tiles is needed.
! Also include effective conductance over non-PFT files for diagnostics.
      IF (sf(462,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(ntiles,len_stlist,                           &
     &    stlist(1,stindex(1,462,3,im_index)),                          &
     &    plltile,stash_pseudo_levels,num_stash_pseudo,                  &
     &    icode,cmessage)
        IF (icode >  0) THEN
          cmessage=                                                     &
     &    "imp_ctl  : error in set_pseudo_list (item 462 = gc)"
          goto 9999
        END IF
       pslevel_out=0
        DO pslevel=1,ntiles
          IF (plltile(pslevel)) THEN
            pslevel_out = pslevel_out+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS(                                      &
     &        stashwork(si(462,3,im_index)+(pslevel_out-1)              &
     &        *row_length*rows),gc(1,pslevel_out),                      &
     &        land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF


! Canopy water on tiles (Kg/m2).
! ----------------------------------------------------------------------
! Item 321 Tiled canopy water

      IF (SF(321,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,321,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 321 = canopy)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(321,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),canopy(1,PSLEVEL_OUT),               &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Canopy capacity on tiles (Kg/m2).
! ----------------------------------------------------------------------
! Item 322 Tiled canopy capacity

      IF (SF(322,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,322,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 322 = catch)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(322,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),catch(1,PSLEVEL_OUT),                &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Roughness length on tiles (M).
! ----------------------------------------------------------------------
! Item 324 Tiled z0

      IF (SF(324,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,324,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 324 = z0m_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(324,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),z0m_tile(1,PSLEVEL_OUT),             &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! 1.5m temperature over tiles
! ----------------------------------------------------------------------
! Item 328 Tiled T at 1.5m

      IF (SF(328,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,328,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 328 = t1p5m_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(328,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),t1p5m_tile(1,PSLEVEL_OUT),           &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! 1.5m specific humidity over tiles
! ----------------------------------------------------------------------
! Item 329 Tiled Q at 1.5 m

      IF (SF(329,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,329,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 329 = q1p5m_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(329,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),q1p5m_tile(1,PSLEVEL_OUT),           &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Surface latent heat flux on tiles (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 330 Tiled latent heat flux

      IF (SF(330,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,330,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 330 = le_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(330,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),le_tile(1,PSLEVEL_OUT),              &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Sublimation on tiles
! ----------------------------------------------------------------------
! Item 331 Tiled sublimation

      IF (SF(331,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,331,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 331 = ei_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(331,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),ei_tile(1,PSLEVEL_OUT),              &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Gross primary productiviy on plant functional types (Kg C/m2/s)
! ----------------------------------------------------------------------
! Item 289 GPP on PFTs
! kdcorbin, 11/10 - changed to tiles 

      IF (SF(289,3)) THEN

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,289,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 289 = gpp_ft)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(289,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),gpp_ft(1,PSLEVEL_OUT),               &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Net primary productiviy on plant functional types (Kg C/m2/s)
! ----------------------------------------------------------------------
! Item 291 NPP on PFTs
! kdcorbin, 11/10 - changed from PFTs to Tiles

      IF (SF(291,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,291,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 291 = npp_ft)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(291,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),npp_ft(1,PSLEVEL_OUT),               &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Plant respiration on plant functional types (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 292 Plant respiration on PFTs
! kdcorbin, 11/10 - changed from PFTs to Tiles

      IF (SF(292,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,292,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 292 = resp_p_ft)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(292,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),resp_p_ft(1,PSLEVEL_OUT),            &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Soil moisture availability factor on plant functional types
! ----------------------------------------------------------------------
! Item 313 FSMC on PFTs

      IF (SF(313,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,313,3,im_index)),                       &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 313 = fsmc)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(313,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),fsmc(1,PSLEVEL_OUT),                 &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF



! Leaf turnover rate on plant functional types
! ----------------------------------------------------------------------
! Item 325 Leaf turnover rate on PFTs
! kdcorbin, 11/10 - changed from PFTs to tiles

      IF (SF(325,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,325,3,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 325 = g_leaf)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1

! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(325,3,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),g_leaf(1,PSLEVEL_OUT),               &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF


! CO2 total flux to atmosphere
! ----------------------------------------------------------------------
! Item 327 CO2 flux

        If (sf(327,3)) Then

! DEPENDS ON: copydiag
           Call copydiag(STASHwork(si(327,3,im_index)),                 &
     &          co2_flux_tot,                                           &
     &          row_length,rows,0,0,0,0, at_extremity,                  &
     &          atmos_im,3,327,                                         &
     &          icode,cmessage)

           If (icode  >   0) Then
              cmessage=                                                 &
     &         "imp_ctl  : error in copydiag(item 327 = co2_flux_tot)"
              goto 9999
           Endif

        End if


! CO2 land surface flux
! ----------------------------------------------------------------------
! Item 326 Land CO2

        IF (SF(326,3)) THEN
! DEPENDS ON: from_land_points
          CALL FROM_LAND_POINTS (                                       &
     &        STASHWORK(SI(326,3,im_index)),                            &
     &         land_co2,                                                &
     &         land_sea_mask,row_length*rows,land_points_dum)
        END IF

! Lestevens
! Net Surface Radiation (W/m2)
! ----------------------------------------------------------------------
! Item 333

      If (sf(333,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(333,3,im_index)),                   &
     &        surf_radflux,                                             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,333,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=                                                   &
     &       "imp_ctl  : error in copydiag(item 333=surf_radflux)"
            goto 9999
         Endif

      End if

! Land heat flux from surface to level 1 (land mean) (W/m2)
! ----------------------------------------------------------------------
! Item 337

      If (sf(337,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(337,3,im_index)),                   &
     &        surf_ht_flux_land,                                        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,337,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=                                                   &
     &       "imp_ctl  : error in copydiag(item 337=surf_ht_flux_land)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
! Net surface sea-ice heat flux (sea mean) (W/m2)
! ----------------------------------------------------------------------
! Item 338

      If (sf(338,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(338,3,im_index)),                   &
     &        surf_ht_flux_sice,                                        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,338,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=                                                   &
     &       "imp_ctl  : error in copydiag(item 338=surf_ht_flux_sice)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!   RIB - measure of surface stability
! ----------------------------------------------------------------------
! Item 339 rib over meaned over open sea and sea-ice

      If (sf(339,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(339,3,im_index)),                   &
     &        rib_ssi,                                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,339,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(rib_ssi)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
! Land Mean 1.5m air temperature
! ----------------------------------------------------------------------
! Item 341 Land mean 1.5m air temperature

      If (sf(341,3)) Then

         do l=1,land_points
           t1p5m_land(l)=0.0
           do n=1,ntiles
             t1p5m_land(l)=t1p5m_land(l)+tile_frac(l,n)*t1p5m_tile(l,n)
           enddo
         enddo

! DEPENDS ON: from_land_points
         CALL FROM_LAND_POINTS (                                        &
     &      STASHWORK(SI(341,3,im_index)),                              &
     &       t1p5m_land,                                                &
     &       land_sea_mask,row_length*rows,land_points_dum)

      End if

! ----------------------------------------------------------------------
! Land Mean 1.5m Specific Humidity
! ----------------------------------------------------------------------
! Item 342 Land mean 1.5m specific humidity

      If (sf(342,3)) Then

         do l=1,land_points
           q1p5m_land(l)=0.0
           do n=1,ntiles
             q1p5m_land(l)=q1p5m_land(l)+tile_frac(l,n)*q1p5m_tile(l,n)
           enddo
         enddo

! DEPENDS ON: from_land_points
         CALL FROM_LAND_POINTS (                                        &
     &      STASHWORK(SI(342,3,im_index)),                              &
     &       q1p5m_land,                                                &
     &       land_sea_mask,row_length*rows,land_points_dum)

      End if
! Surface sensible heat flux  meaned over open sea and sea-ice
! ----------------------------------------------------------------------
! Item 343
! surface = ftl(,,1)

      If (sf(343,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(343,3,im_index)),                   &
     &        ftl_ssi,                                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,343,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(ftl_ssi)"
            goto 9999
         Endif

      End if
! Surface moisture flux  meaned over open sea and sea-ice
! ----------------------------------------------------------------------
! Item 347
! surface = fqw(,,1)
      If (sf(347,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(347,3,im_index)),                   &
     &        e_ssi,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,347,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(e_ssi)"
            goto 9999
         Endif

      End if


! Sublimation over sea-ice meaned over sea portion of gridbox
! ----------------------------------------------------------------------
! Item 353

      If (sf(353,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(353,3,im_index)),                   &
     &        ei_sice,                                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,353,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(ei_sice)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
! Sea-ice surface net radiation (sea mean) (W/m2)
! ----------------------------------------------------------------------
! Item 381

      If (sf(381,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(381,3,im_index)),                   &
     &        radnet_sice,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,381,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=                                                   &
     &       "imp_ctl  : error in copydiag(item 381=radnet_sice)"
            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
! Land fraction
! ----------------------------------------------------------------------
! Item 395

      If (sf(395,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(395,3,im_index)),                   &
     &        flandg,                                                   &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,395,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=                                                   &
     &       "imp_ctl  : error in copydiag(item 395=flandg)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
! Wind shear over land
! ----------------------------------------------------------------------
! Item 389

      If (sf(389,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(389,3,im_index)),                   &
     &        vshr_land,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,389,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(vshr_land)"
            goto 9999
         Endif

      End if

! Wind shear over sea
! ----------------------------------------------------------------------
! Item 390

      If (sf(390,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(390,3,im_index)),                   &
     &        vshr_ssi,                                                 &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,390,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(vshr_ssi)"
            goto 9999
         Endif

      End if

! X-component of surface stress
! ----------------------------------------------------------------------
! Item 460

      If (sf(460,3)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(460,3,im_index)),                   &
     &        taux(1:row_length, 1:rows, 1),                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,460,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(taux_surf)"
            goto 9999
         Endif

      End if

! Y-component of surface stress
! ----------------------------------------------------------------------
! Item 461

      If (sf(461,3)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(461,3,im_index)),                   &
     &        tauy(1:row_length, 1:n_rows, 1),                          &
     &        row_length,n_rows,0,0,0,0, at_extremity,                  &
     &        atmos_im,3,461,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(tauy_surf)"
            goto 9999
         Endif

      End if

! X-component of stress over land
! ----------------------------------------------------------------------
! Item 391

      If (sf(391,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(391,3,im_index)),                   &
     &        taux_land,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,391,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(taux_land)"
            goto 9999
         Endif

      End if

! X-component of stress over sea
! ----------------------------------------------------------------------
! Item 392

      If (sf(392,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(392,3,im_index)),                   &
     &        taux_ssi,                                                 &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,3,392,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(taux_ssi)"
            goto 9999
         Endif

      End if

! Y-component of stress over land
! ----------------------------------------------------------------------
! Item 393

      If (sf(393,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(393,3,im_index)),                   &
     &        tauy_land,                                                &
     &        row_length,n_rows,0,0,0,0, at_extremity,                  &
     &        atmos_im,3,393,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(tauy_land)"
            goto 9999
         Endif

      End if

! Y-component of stress over sea
! ----------------------------------------------------------------------
! Item 394

      If (sf(394,3)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(394,3,im_index)),                   &
     &        tauy_ssi,                                                 &
     &        row_length,n_rows,0,0,0,0, at_extremity,                  &
     &        atmos_im,3,394,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(tauy_ssi)"
            goto 9999
         Endif

      End if

! Diagnostics required for soil moisture nudging scheme macro

! Item 50, RHOKH, exchange coefficients for moisture

      IF( sf(50,3)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(50,3,im_index)),rhokh,           &
     &        row_length,rows,boundary_layer_levels,0,0,0,0,            &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,50,3,im_index)),len_stlist,            &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,3,50,                                            &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="imp_ctl  : error in copydiag(RHOKH)"
            goto 9999
         ENDIF

      ENDIF
!
! Item 51, RESFS, Combined soil, stomatol and aerodynamic
!                 resistance to evaporation
      IF (SF(51,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,51,3,im_index)),                        &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 51 = resfs)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(51,3,im_index)+(PSLEVEL_OUT-1)             &
     &           *row_length*rows),resfs(1,PSLEVEL_OUT),                &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF
!
! Item 52, CHR1P5M, Ratio of coefficients required for calculation
!                 of 1.5m temperature
      IF (SF(52,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,52,3,im_index)),                        &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 52 = chr1p5m)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(52,3,im_index)+(PSLEVEL_OUT-1)             &
     &           *row_length*rows),chr1p5m(1,PSLEVEL_OUT),              &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF

! Item 53, ALPHA1, Gradient of saturated humidity with respect to
!                  tempearture between surface and model level 1
!                 of 1.5m temperature
      IF (SF(53,3)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,53,3,im_index)),                        &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "imp_ctl  : error in set_pseudo_list(item 53 = alpha1)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(53,3,im_index)+(PSLEVEL_OUT-1)             &
     &           *row_length*rows),alpha1(1,PSLEVEL_OUT),               &
     &           land_sea_mask,row_length*rows,land_points_dum)
          END IF
        END DO
      END IF

! Item 54, RA Aerodyanmic resiatance (s/m)
! ----------------------------------------------------------------------

      IF (SF(54,3)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(54,3,im_index)),                               &
     &       ra,                                                        &
     &       land_sea_mask,row_length*rows,land_points_dum)
      END IF


! Item 55, WT_EXT
! Culmulative transpiration from soil (needed for soil moisture nudging)
! ----------------------------------------------------------------------

      If (sf(55,3)) Then

         Do k = 1, dsm_levels
            Do j= 1, rows
               Do i = 1, row_length
                 l = i + (j-1)*row_length + (k-1)*rows*row_length
                 interp_data_3(l) = rmdi
               End Do
            End Do
            Do i = 1, land_points
               l = land_index(i) + (k-1)*rows*row_length
               interp_data_3(l) = wt_ext(i,k)
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(55,3,im_index)),interp_data_3,   &
     &        row_length,rows,dsm_levels,0,0,0,0,                       &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,55,3,im_index)),len_stlist,            &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,3,55,                                            &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(item 55)"
         Endif

      End if
!
! MGS. Extra boundary layer diagnostics required for UKCA model
! -------------------------------------------------------------------
!
! Item 60, RHOKH_MIX

      IF( sf(60,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(60,sect,im_index)),              &
     &        rhokh_mix,                                                &
     &        row_length,rows,boundary_layer_levels,0,0,0,0,            &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,60,sect,im_index)),len_stlist,         &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,60,                                         &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="imp_ctl  : error in copydiag(RHOKH_MIX)"
            goto 9999
         ENDIF

      ENDIF
!
! Item 61, RHO_ARESIST

      If (sf(61,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(61,sect,im_index)),                 &
     &        rho_aresist,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,61,                                         &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(rho_aresist)"
            goto 9999
         Endif

      End if
!
! Item 62, ARESIST [ 1/(CD_STD*VSHR) ]

      If (sf(62,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(62,sect,im_index)),                 &
     &        aresist,                                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,62,                                         &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode > 0) Then
            cmessage="imp_ctl  : error in copydiag(aresist)"
            goto 9999
         Endif

      End if
!
! Item 63, RESIST_B

      If (sf(63,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(63,sect,im_index)),                 &
     &        resist_b,                                                 &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,63,                                         &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(resist_b)"
            goto 9999
         Endif

      End if
!
! Item 64, DTRDZ_CHARNEY_GRID

      IF( sf(64,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(64,sect,im_index)),              &
     &        dtrdz_charney_grid,                                       &
     &        row_length,rows,boundary_layer_levels,0,0,0,0,            &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,64,sect,im_index)),len_stlist,         &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,64,                                         &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="imp_ctl  : error in copydiag(dtrdz_charney_grid)"
            goto 9999
         ENDIF

      ENDIF
!
! Item 65, GRID-LEVEL OF SML INVERSION

      If (sf(65,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(65,sect,im_index)),                 &
     &        kent,                                                     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,65,                                         &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="imp_ctl  : error in copydiag(tauy_ssi)"
            goto 9999
         Endif

      End if

! IF any of diagnsotics 66-68, 70-72 are selected, set interp_data_bl
! to zero
      IF (sf(66,sect) .OR. sf(67,sect) .OR. sf(68,sect) .OR.            &
     &    sf(70,sect) .OR. sf(71,sect) .OR. sf(72,sect)) THEN
        DO k = 1, npft
          DO j = 1, rows
            DO i = 1, row_length
              interp_data_bl(i,j,k) = 0.0
            END DO
          END DO
        END DO

      END IF
!
! Item 66, Rho * entrainment rate  (we_lim)

      IF (sf(66,sect)) THEN

        DO k = 1, 3
          DO j = 1, rows
            DO i = 1, row_length
              interp_data_bl(i,j,k) = we_lim(i,j,k)
            END DO
          END DO
        END DO

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(npft,len_stlist,                           &
     &       stlist(1,stindex(1,66,sect,im_index)),                     &
     &       pllbl,stash_pseudo_levels,num_stash_pseudo,                &
     &       icode,cmessage)
        IF (icode >  0) THEN
          cmessage =                                                    &
     &    "imp_ctl  : error in set_pseudo_list(item 66 = we_lim)"
          GOTO 9999
        END IF
        pslevel_out = 0
        DO pslevel = 1, npft
          IF (pllbl(pslevel)) THEN
            pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
            Call copydiag(STASHwork(si(66,sect,im_index)+(pslevel_out-1)&
     &        *row_length*rows),                                        &
     &        interp_data_bl(1,1,pslevel),                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,66,                                         &
     &        icode,cmessage)
            If (icode  >   0) Then
              cmessage=": error in copydiag(item 66)"
            Endif
          END IF
        END DO
      END IF
!
! Item 67, Fraction of the timestep (t_frac)

      IF (sf(67,sect)) THEN

        DO k = 1, 3
          DO j = 1, rows
            DO i = 1, row_length
              interp_data_bl(i,j,k) = t_frac(i,j,k)
            END DO
          END DO
        END DO

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(npft,len_stlist,                           &
     &       stlist(1,stindex(1,67,sect,im_index)),                     &
     &       pllbl,stash_pseudo_levels,num_stash_pseudo,                &
     &       icode,cmessage)
        IF (icode >  0) THEN
          cmessage =                                                    &
     &    "imp_ctl  : error in set_pseudo_list(item 67 = t_frac)"
          GOTO 9999
        END IF

        pslevel_out = 0
        DO pslevel = 1, npft
          IF (pllbl(pslevel)) THEN
            pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
            Call copydiag(STASHwork(si(67,sect,im_index)+(pslevel_out-1)&
     &        *row_length*rows),                                        &
     &        interp_data_bl(1,1,pslevel),                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,67,                                         &
     &        icode,cmessage)
            If (icode  >   0) Then
              cmessage=": error in copydiag(item 67)"
            Endif
          END IF
        END DO
      END IF
!
! Item 68, zrzi

      IF( sf(68,sect)) THEN

        DO k = 1, 3
          DO j = 1, rows
            DO i = 1, row_length
              interp_data_bl(i,j,k) = zrzi(i,j,k)
            END DO
          END DO
        END DO

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(npft,len_stlist,                           &
     &       stlist(1,stindex(1,68,sect,im_index)),                     &
     &       pllbl,stash_pseudo_levels,num_stash_pseudo,                &
     &       icode,cmessage)
        IF (icode >  0) THEN
          cmessage =                                                    &
     &    "imp_ctl  : error in set_pseudo_list(item 68 = zrzi)"
          GOTO 9999
        END IF

        pslevel_out = 0
        DO pslevel = 1, npft
          IF (pllbl(pslevel)) THEN
            pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
            Call copydiag(STASHwork(si(68,sect,im_index)+(pslevel_out-1)&
     &        *row_length*rows),                                        &
     &        interp_data_bl(1,1,pslevel),                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,68,                                         &
     &        icode,cmessage)
            If (icode  >   0) Then
              cmessage=": error in copydiag(item 68)"
            Endif
          END IF
        END DO
      ENDIF
!
! Item 69, GRID-LEVEL OF DSC INVERSION (kent_dsc)

      If (sf(69,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(69,sect,im_index)),                 &
     &        kent_dsc,                                                 &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,69,                                         &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode > 0) Then
            cmessage="imp_ctl  : error in copydiag(kent_dsc)"
            goto 9999
         Endif

      End if
!
! Item 70, Rho * entrainment rate  (we_lim_dsc)

      IF( sf(70,sect)) THEN

        DO k = 1, 3
          DO j = 1, rows
            DO i = 1, row_length
              interp_data_bl(i,j,k) = we_lim_dsc(i,j,k)
            END DO
          END DO
        END DO

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(npft,len_stlist,                           &
     &       stlist(1,stindex(1,70,sect,im_index)),                     &
     &       pllbl,stash_pseudo_levels,num_stash_pseudo,                &
     &       icode,cmessage)
        IF (icode >  0) THEN
          cmessage =                                                    &
     &    "imp_ctl  : error in set_pseudo_list(item 70 = we_lim_dsc)"
          GOTO 9999
        END IF

        pslevel_out = 0
        DO pslevel = 1, npft
          IF (pllbl(pslevel)) THEN
            pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
            Call copydiag(STASHwork(si(70,sect,im_index)+(pslevel_out-1)&
     &        *row_length*rows),                                        &
     &        interp_data_bl(1,1,pslevel),                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,70,                                         &
     &        icode,cmessage)
            If (icode  >   0) Then
              cmessage=": error in copydiag(item 70)"
            Endif
          END IF
        END DO
      ENDIF
!
! Item 71, Fraction of the timestep (t_frac_dsc)

      IF (sf(71,sect)) THEN

        DO k = 1, 3
          DO j = 1, rows
            DO i = 1, row_length
              interp_data_bl(i,j,k) = t_frac_dsc(i,j,k)
            END DO
          END DO
        END DO

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(npft,len_stlist,                           &
     &       stlist(1,stindex(1,71,sect,im_index)),                     &
     &       pllbl,stash_pseudo_levels,num_stash_pseudo,                &
     &       icode,cmessage)
        IF (icode >  0) THEN
          cmessage =                                                    &
     &    "imp_ctl  : error in set_pseudo_list(item 71 = t_frac_dsc)"
          GOTO 9999
        END IF

        pslevel_out = 0
        DO pslevel = 1, npft
          IF (pllbl(pslevel)) THEN
            pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
            Call copydiag(STASHwork(si(71,sect,im_index)+(pslevel_out-1)&
     &        *row_length*rows),                                        &
     &        interp_data_bl(1,1,pslevel),                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,71,                                         &
     &        icode,cmessage)
            If (icode > 0) Then
              cmessage=": error in copydiag(item 71)"
            Endif
          END IF
        END DO
      ENDIF
!
! Item 72, zrzi_dsc

      IF (sf(72,sect)) THEN

        DO k = 1, 3
          DO j = 1, rows
            DO i = 1, row_length
              interp_data_bl(i,j,k) = zrzi_dsc(i,j,k)
            END DO
          END DO
        END DO

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(npft,len_stlist,                           &
     &       stlist(1,stindex(1,72,sect,im_index)),                     &
     &       pllbl,stash_pseudo_levels,num_stash_pseudo,                &
     &       icode,cmessage)
        IF (icode >  0) THEN
          cmessage =                                                    &
     &    "imp_ctl  : error in set_pseudo_list(item 72 = zrzi_dsc)"
          GOTO 9999
        END IF

        pslevel_out = 0
        DO pslevel = 1, npft
          IF (pllbl(pslevel)) THEN
            pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
            Call copydiag(STASHwork(si(72,sect,im_index)+(pslevel_out-1)&
     &        *row_length*rows),                                        &
     &        interp_data_bl(1,1,pslevel),                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,72,                                         &
     &        icode,cmessage)
            If (icode > 0) Then
              cmessage=": error in copydiag(item 72)"
            Endif
          END IF
        END DO
      ENDIF
!
! Item 73, ZHSC  Height of top of decoupled layer

      If (sf(73,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(73,sect,im_index)),                 &
     &        zhsc,                                                     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,73,                                         &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode > 0) Then
            cmessage="imp_ctl  : error in copydiag(zhsc)"
            goto 9999
         Endif

      End if
!
! Items 74-79 Surface layer resist for dust

      Do k = 1, ndiv
        item = 73 + k
        If (sf(item,sect)) Then
          Do j = 1, rows
            Do i = 1, row_length
              interp_data(i,j) = r_b_dust(i,j,k)
            End Do
          End Do

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),              &
     &        interp_data,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
!#include <argppx/argppx.h>
     &        icode,cmessage)

          If (icode > 0) Then
            cmessage="imp_ctl  : error in copydiag(r_b_dust div  )"
            write(cmessage(43:43),'(I1)') k
            goto 9999
          Endif

        End if
      End do

! Mineral Dust Diagnostics
! -------------------------------------------------------------------
!
      IF (L_DUST) THEN

        DO I = 1,NDIV
          ITEM = 400+I  ! dust emission flux
          IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN

! DEPENDS ON: copydiag
            CALL COPYDIAG (STASHWORK(SI(ITEM,SECT,IM_INDEX)),           &
     &        DUST_FLUX(1:ROW_LENGTH,1:ROWS,I),                         &
     &        ROW_LENGTH,ROWS,0,0,0,0, AT_EXTREMITY,                    &
     &        ATMOS_IM,SECT,ITEM,                                       &
     &        ICODE,CMESSAGE)
            IF (ICODE >  0) THEN
              CMESSAGE=": ERROR IN COPYDIAG_3D(ITEM 401-6)"//CMESSAGE
            ENDIF

          ENDIF  !  SF(ITEM,SECT)
        ENDDO !NDIV
!
!
!
        DO I = 1,NDIV
          ITEM = 410+I !threshold friction velocity on tiles
          IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN
! DEPENDS ON: set_pseudo_list
            CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                     &
     &       STLIST(1,STINDEX(1,ITEM,SECT,IM_INDEX)),                   &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
            IF (ICODE >  0) THEN
              CMESSAGE=                                                 &
     &        "IMP_CTL  : ERROR IN SET_PSEUDO_LIST(ITEM 411-416)"
            ENDIF

            PSLEVEL_OUT=0
            DO PSLEVEL=1,NTILES
              IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                 &
     &           STASHWORK(SI(ITEM,SECT,IM_INDEX)+(PSLEVEL_OUT-1)       &
     &           *ROW_LENGTH*ROWS),U_S_T_TILE(1,PSLEVEL_OUT,I),         &
     &           LAND_SEA_MASK,ROW_LENGTH*ROWS,LAND_POINTS_DUM)
              ENDIF
            ENDDO !NTILES
          ENDIF
        ENDDO !NDIV
!
!
!
        DO I = 1,NDIV
          ITEM = 420+I !dry threshold friction velocity on tiles
          IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN
! DEPENDS ON: set_pseudo_list
            CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                     &
     &       STLIST(1,STINDEX(1,ITEM,SECT,IM_INDEX)),                   &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
            IF (ICODE >  0) THEN
              CMESSAGE=                                                 &
     &        "IMP_CTL  : ERROR IN SET_PSEUDO_LIST(ITEM 421-426)"
            ENDIF

            PSLEVEL_OUT=0
            DO PSLEVEL=1,NTILES
              IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                 &
     &           STASHWORK(SI(ITEM,SECT,IM_INDEX)+(PSLEVEL_OUT-1)       &
     &           *ROW_LENGTH*ROWS),U_S_T_DRY_TILE(1,PSLEVEL_OUT,I),     &
     &           LAND_SEA_MASK,ROW_LENGTH*ROWS,LAND_POINTS_DUM)
               ENDIF
             END DO !NTILES
           ENDIF
         ENDDO !NDIV

!
!
!
         ITEM = 430 ! friction velocity on tiles
         IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN
! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                      &
     &       STLIST(1,STINDEX(1,ITEM,SECT,IM_INDEX)),                   &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             CMESSAGE=                                                  &
     &       "IMP_CTL  : ERROR IN SET_PSEUDO_LIST(ITEM 431-436)"
           ENDIF

           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
               PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
               CALL FROM_LAND_POINTS (                                  &
     &          STASHWORK(SI(ITEM,SECT,IM_INDEX)+(PSLEVEL_OUT-1)        &
     &          *ROW_LENGTH*ROWS),U_S_STD_TILE(1,PSLEVEL_OUT),          &
     &          LAND_SEA_MASK,ROW_LENGTH*ROWS,LAND_POINTS_DUM)
             ENDIF
           ENDDO
         ENDIF
!
!
!
         DO I = 1,NDIV
          ITEM = 450+I  ! dry deposition from 2nd layer
          IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN

! DEPENDS ON: copydiag
            CALL COPYDIAG (STASHWORK(SI(ITEM,SECT,IM_INDEX)),           &
     &        DRYDEP2(1:ROW_LENGTH,1:ROWS,I),                           &
     &        ROW_LENGTH,ROWS,0,0,0,0, AT_EXTREMITY,                    &
     &        ATMOS_IM,SECT,ITEM,                                       &
     &        ICODE,CMESSAGE)
            IF (ICODE >  0) THEN
              CMESSAGE=": ERROR IN COPYDIAG_3D(ITEM 451-6)"//CMESSAGE
            ENDIF

          ENDIF  !  SF(ITEM,SECT)
        ENDDO !NDIV

      ENDIF !L_DUST

! CABLE case
! Lest May 2011 - changed vcodes from 140s to 800s
! Soil temperature on tiles
! ----------------------------------------------------------------------

! Items 801-804/801-806 Tiled TSOIL layers 1-4/1-6

      do k=1,dsm_levels
        vcode=800+k
        IF (SF(vcode,3)) THEN
           ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item ", vcode, " = tsoil_tile)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
      &           *row_length*rows),tsoil_tile(1,PSLEVEL_OUT,k),           &
      &           land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
         END IF
       end do

! Soil moisture on tiles
! ----------------------------------------------------------------------
! Les 2010 - changed vcodes to allow for 6 layers
! Items 805-808/807-812 Tiled SMCL layers 1-4/1-6

      do k=1,dsm_levels
!        vcode=144+k
        vcode=800+dsm_levels+k
        IF (SF(vcode,3)) THEN
           ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item ", vcode, " = smcl_tile)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
      &           *row_length*rows),smcl_tile(1,PSLEVEL_OUT,k),           &
      &           land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
         END IF
       end do

! Soil moisture frozen fraction on tiles
! ----------------------------------------------------------------------

! Items 809-812/813-818 Tiled STHF layers 1-4/1-6

      do k=1,dsm_levels
!        vcode=148+k
        vcode=800+(2*dsm_levels)+k
        IF (SF(vcode,3)) THEN
           ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item ", vcode, " = sthf_tile)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
      &           *row_length*rows),sthf_tile(1,PSLEVEL_OUT,k),           &
      &           land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
         END IF
       end do

! Snow depth on tiles
! ----------------------------------------------------------------------

! Items 813-815/819-821(if 6 soil layers) Tiled SNOW_DEPTH3L layers 1-3

      do k=1,3
!        vcode=152+k
        vcode=800+(3*dsm_levels)+k
        IF (SF(vcode,3)) THEN
           ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item ", vcode, " = snowdepth_tile)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
      &           *row_length*rows),snow_depth3l(1,PSLEVEL_OUT,k),         &
      &           land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
         END IF
       end do

! Snow mass on tiles
! ----------------------------------------------------------------------

! Items 816-818/822-824(if 6 soil layers) Tiled SNOW_MASS3L layers 1-3

      do k=1,3
!        vcode=155+k
        vcode=800+(3*dsm_levels+3)+k
        IF (SF(vcode,3)) THEN
           ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item ", vcode, " = snowmass_tile)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
      &           *row_length*rows),snow_mass3l(1,PSLEVEL_OUT,k),         &
      &           land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
         END IF
       end do

! Snow temperature on tiles
! ----------------------------------------------------------------------

! Items 819-821/825-827(if 6 soil layers) Tiled SNOW_TMP3L layers 1-3

      do k=1,3
!        vcode=158+k
        vcode=800+(3*dsm_levels+6)+k
        IF (SF(vcode,3)) THEN
           ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item ", vcode, " = snowtmp_tile)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
      &           *row_length*rows),snow_tmp3l(1,PSLEVEL_OUT,k),           &
      &           land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
         END IF
       end do

! Snow density on tiles
! ----------------------------------------------------------------------

! Items 822-824/828-830(if 6 soil layers) Tiled SNOW_RHO3L layers 1-3

      do k=1,3
!        vcode=161+k
        vcode=800+(3*dsm_levels+9)+k
        IF (SF(vcode,3)) THEN
           ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item ", vcode, " = snowrho3l_tile)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
      &           *row_length*rows),snow_rho3l(1,PSLEVEL_OUT,k),           &
      &           land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
         END IF
       end do

! Mean snow density
! ----------------------------------------------------------------------

! Items 825/831(if 6 soil layers) Tiled SNOW_RHO1L 

!      vcode=165
      vcode=800+(3*dsm_levels+12)+1
      IF (SF(vcode,3)) THEN
         ! DEPENDS ON: set_pseudo_list
         CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &      STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &      PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &      ICODE,CMESSAGE)
         IF (ICODE >  0) THEN
           write(cmessage,"(a,i3.3,a)")                                  &
     &       "imp_ctl  : error in set_pseudo_list(item ", vcode, " = snow_rho1l)"
           goto 9999
         END IF
         PSLEVEL_OUT=0
         DO PSLEVEL=1,NTILES
           IF (PLLTILE(PSLEVEL)) THEN
              PSLEVEL_OUT=PSLEVEL_OUT+1
              ! DEPENDS ON: from_land_points
              CALL FROM_LAND_POINTS (                                    &
              STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
                *row_length*rows),snow_rho1l(1,PSLEVEL_OUT),             &
                land_sea_mask,row_length*rows,land_points_dum)
             END IF
          END DO
       END IF

 ! Snow age
! ----------------------------------------------------------------------

! Items 826/832(if 6 soil layers) Tiled SNAGE_TILE

!      vcode=166
      vcode=800+(3*dsm_levels+13)+1
      IF (SF(vcode,3)) THEN
         ! DEPENDS ON: set_pseudo_list
         CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &      STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &      PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &      ICODE,CMESSAGE)
         IF (ICODE >  0) THEN
           write(cmessage,"(a,i3.3,a)")                                  &
     &       "imp_ctl  : error in set_pseudo_list(item ", vcode, " = snage_tile)"
           goto 9999
         END IF
         PSLEVEL_OUT=0
         DO PSLEVEL=1,NTILES
           IF (PLLTILE(PSLEVEL)) THEN
              PSLEVEL_OUT=PSLEVEL_OUT+1
              ! DEPENDS ON: from_land_points
              CALL FROM_LAND_POINTS (                                    &
              STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
                *row_length*rows),snage_tile(1,PSLEVEL_OUT),             &
                land_sea_mask,row_length*rows,land_points_dum)
             END IF
          END DO
       END IF

! CASA-CNP case
! Lest Sept 2012 - add casaCNP diagnostics
! Carbon Pools on Tiles
! ----------------------------------------------------------------------

! Items 851-860 Tiled CPOOL_TILE

      do k=1,10
        vcode=850+k
        IF (SF(vcode,3)) THEN
           ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item ", vcode, " = cpool_tile)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
      &           *row_length*rows),cpool_tile(1,PSLEVEL_OUT,k),           &
      &           land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
         END IF
       end do

! Nitrogen Pools on Tiles
! ----------------------------------------------------------------------

! Items 861-870 Tiled NPOOL_TILE

      do k=1,10
        vcode=860+k
        IF (SF(vcode,3)) THEN
           ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item ", vcode, " = npool_tile)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
      &           *row_length*rows),npool_tile(1,PSLEVEL_OUT,k),           &
      &           land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
         END IF
       end do

! Phosphorus Pools on Tiles
! ----------------------------------------------------------------------

! Items 871-882 Tiled PPOOL_TILE

      do k=1,12
        vcode=870+k
        IF (SF(vcode,3)) THEN
           ! DEPENDS ON: set_pseudo_list
           CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &        STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &        PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &        ICODE,CMESSAGE)
           IF (ICODE >  0) THEN
             write(cmessage,"(a,i3.3,a)")                                  &
     &         "imp_ctl  : error in set_pseudo_list(item ", vcode, " = ppool_tile)"
             goto 9999
           END IF
           PSLEVEL_OUT=0
           DO PSLEVEL=1,NTILES
             IF (PLLTILE(PSLEVEL)) THEN
                PSLEVEL_OUT=PSLEVEL_OUT+1
                ! DEPENDS ON: from_land_points
                CALL FROM_LAND_POINTS (                                    &
      &         STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
      &           *row_length*rows),ppool_tile(1,PSLEVEL_OUT,k),           &
      &           land_sea_mask,row_length*rows,land_points_dum)
               END IF
            END DO
         END IF
       end do

! Soil Order
! ----------------------------------------------------------------------

! Item 883 SOIL_ORDER

      If (sf(883,3)) Then

     !    do l=1,land_points
     !      t1p5m_land(l)=0.0
     !      do n=1,ntiles
     !        t1p5m_land(l)=t1p5m_land(l)+tile_frac(l,n)*t1p5m_tile(l,n)
     !      enddo
     !    enddo

! DEPENDS ON: from_land_points
         CALL FROM_LAND_POINTS (                                        &
     &      STASHWORK(SI(883,3,im_index)),                              &
     &       soil_order,                                                &
     &       land_sea_mask,row_length*rows,land_points_dum)

      End if

! Transpiration
! ----------------------------------------------------------------------
! Les Sept 2012
! Items 174 Tiled TRANSP_TILE

!     vcode=166
      vcode=140+(3*dsm_levels+15)+1
      IF (SF(vcode,3)) THEN
         ! DEPENDS ON: set_pseudo_list
         CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &      STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &      PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &      ICODE,CMESSAGE)
         IF (ICODE >  0) THEN
           write(cmessage,"(a,i3.3,a)")                                  &
     &       "imp_ctl  : error in set_pseudo_list(item ", vcode, " = transp_tile)"
           goto 9999
         END IF
         PSLEVEL_OUT=0
         DO PSLEVEL=1,NTILES
           IF (PLLTILE(PSLEVEL)) THEN
              PSLEVEL_OUT=PSLEVEL_OUT+1
              ! DEPENDS ON: from_land_points
              CALL FROM_LAND_POINTS (                                    &
              STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
                *row_length*rows),transp_tile(1,PSLEVEL_OUT),            &
                land_sea_mask,row_length*rows,land_points_dum)
             END IF
          END DO
       END IF

! Leaf Area Index from Casa-CNP (GLAI)
! ----------------------------------------------------------------------

! Items 893 (if 6 soil layers) Tiled GLAI

      vcode=893
!      vcode=140+(3*dsm_levels+15)+1
      IF (SF(vcode,3)) THEN
         ! DEPENDS ON: set_pseudo_list
         CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &      STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &      PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &      ICODE,CMESSAGE)
         IF (ICODE >  0) THEN
           write(cmessage,"(a,i3.3,a)")                                  &
     &       "imp_ctl  : error in set_pseudo_list(item ", vcode, " = glai)"
           goto 9999
         END IF
         PSLEVEL_OUT=0
         DO PSLEVEL=1,NTILES
           IF (PLLTILE(PSLEVEL)) THEN
              PSLEVEL_OUT=PSLEVEL_OUT+1
              ! DEPENDS ON: from_land_points
              CALL FROM_LAND_POINTS (                                    &
              STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
                *row_length*rows),glai(1,PSLEVEL_OUT),            &
                land_sea_mask,row_length*rows,land_points_dum)
             END IF
          END DO
       END IF
!
! Phenology Phase from Casa-CNP (PHENPHASE)
! ----------------------------------------------------------------------

! Items 894 (if 6 soil layers) Tiled PHENPHASE

      vcode=894
!      vcode=140+(3*dsm_levels+15)+1
      IF (SF(vcode,3)) THEN
         ! DEPENDS ON: set_pseudo_list
         CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &      STLIST(1,STINDEX(1,vcode,3,im_index)),                       &
     &      PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,                &
     &      ICODE,CMESSAGE)
         IF (ICODE >  0) THEN
           write(cmessage,"(a,i3.3,a)")                                  &
     &       "imp_ctl  : error in set_pseudo_list(item ", vcode, " = phenphase)"
           goto 9999
         END IF
         PSLEVEL_OUT=0
         DO PSLEVEL=1,NTILES
           IF (PLLTILE(PSLEVEL)) THEN
              PSLEVEL_OUT=PSLEVEL_OUT+1
              ! DEPENDS ON: from_land_points
              CALL FROM_LAND_POINTS (                                    &
              STASHWORK(SI(vcode,3,im_index)+(PSLEVEL_OUT-1)             &
                *row_length*rows),phenphase(1,PSLEVEL_OUT),            &
                land_sea_mask,row_length*rows,land_points_dum)
             END IF
          END DO
       END IF
!
! ----------------------------------------------------------------------
!  single point error handling.
 9999 continue ! to be removed
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
      END SUBROUTINE diagnostics_bl

