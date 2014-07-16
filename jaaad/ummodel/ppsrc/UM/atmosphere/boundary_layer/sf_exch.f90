
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!   SUBROUTINE SF_EXCH------------------------------------------------
!!!
!!!  Purpose: Calculate coefficients of turbulent exchange between
!!!           the surface and the lowest atmospheric layer, and
!!!           "explicit" fluxes between the surface and this layer.
!!!
!!!  Suitable for Single Column use.
!!!
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version 2, dated 18/1/90.
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Project task:
!!!
!!!  Documentation: UM Documentation Paper No 24, section P243.
!!!                 See especially sub-section (ix).
!!!
!!!---------------------------------------------------------------------

!Note by YMA.  Changes are made to alter the momentum roughness length z0m 
!and scalar roughness length z0h.  z0m is changed from a constant value to
!variable values dependent on wind speed at 10 m.  z0h is changed to a new
!expression.  Both the formulae for z0m and z0h are regressed from the 
!observations summarised in Fairall et al., 2003, J. Clim., 571-.

! Arguments :-

      SUBROUTINE SF_EXCH (                                              &
     & ROW_LENGTH,ROWS,OFF_X,OFF_Y,HALO_I,HALO_J,                       &
     & LAND_PTS,NTILES,LAND_INDEX,                                      &
     & TILE_INDEX,TILE_PTS,FLAND,FLANDG,                                &
     & BQ_1,BT_1,CANHC_TILE,CANOPY,CATCH,DZSOIL,FLAKE,GC,HCONS,         &
     & CAN_MODEL,CATCH_SNOW, lq_mix_bl,                                 &
     & HO2R2_OROG,ICE_FRACT,SNOW_TILE,PSTAR,QW_1,RADNET_SICE,           &
     & RADNET_TILE,                                                     &
     & SIL_OROG_LAND,SMVCST,TILE_FRAC,TIMESTEP,                         &
     & TL_1,TI,TS1,COR_UST,COR_MO_ITER,                                 &
     & TSTAR_TILE,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,Z_LAND,     &
     & L_CTILE,SeaSalinityFactor,ISeaZ0T,I_SCRN_T_DIAG,                 &
     & TSTAR,L_spec_z0,Z0M_SCM,Z0H_SCM,L_DUST,                          &
     & VFRAC_TILE,VSHR_LAND,VSHR_SSI,ZH,Z0_TILE,Z1_UV,Z1_TQ,LAND_MASK,  &
     & SU10,SV10,SQ1P5,ST1P5,SFME,SZ0HEFF,LTIMER,FORMDRAG,FD_stab_dep,  &
     & OROG_DRAG_PARAM,Z0MSEA,                                          &
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,CD,CH,                       &
     & RECIP_L_MO_SEA,CDR10M,                                           &
     & CHR1P5M,CHR1P5M_SICE,E_SEA,FME,FQW_1,FQW_TILE,EPOT_TILE,         &
     & FQW_ICE,                                                         &
     & FTL_1,FTL_TILE,FTL_ICE,FRACA,H_BLEND_OROG,H_SEA,Charnock,        &
     & RHOSTAR,RESFS,RESFT,RIB,RIB_TILE,                                &
     & FB_SURF,U_S,Q1_SD,T1_SD,Z0HSSI,Z0H_TILE,Z0H_EFF,                 &
     & Z0M_GB,Z0MSSI,Z0M_TILE,Z0M_EFF,RHO_ARESIST,ARESIST,RESIST_B,     &
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,                     &
     & R_B_DUST,CD_STD_DUST,U_S_STD_TILE,                               &
     & RHO_CD_MODV1,RHOKH_1,RHOKH_1_SICE,RHOKM_1,RHOKM_LAND,RHOKM_SSI,  &
     & RHOKPM,RHOKPM_POT,RHOKPM_SICE,BL_diag,ANTHROP_HEAT               &
! EAK
     &,l_cable                                                          &
     &,SM_LEVELS,N_ROWS,LAND_PTS_TRIF,NPFT_TRIF,DIM_CS1,DIM_CS2         &
     &,SMVCCL,SMVCWT                                                    &
     &,CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE         &
     &,SW_TILE,LW_DOWN,SW_DOWN                                          &
     &,U,V,U_P,V_P                                                      &
     &,surf_down_sw,alb_tile,cos_zenith_angle,l_tile_pts,               &
!     &,surf_down_sw,alb_tile,cos_zenith_angle,               &
     & ls_rain,ls_snow,                                                 &
     & lat,long,day,time_sec,                                           &
!     & )
     & SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                   &
     & SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,             &
     & TSOIL_TILE,T_SURF_TILE,BEXP,                                     &
     & SATHH,SATCON,HCAP,HCON,                                          &
!     & SOIL_TYPE,VEG_TYPE,                                              &
     & SOIL_TYPE,VEG_TYPE,albsoil,                                      &
     & ISNOW_FLG3L,total_nsteps,                                        &
     &            FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,              &
     &            TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,          &
     &            USTAR_CAB,SURF_HTF_CAB,                               &
!     & )
     & SOIL_LAYER_MOISTURE,STHF,STHU,T_SOIL,                            &
     & CS,FRAC,CANHT_FT,PHOTOSYNTH_ACT_RAD,LAI_FT,                      &
     & G_LEAF,GS,                                                       &
     & LE_TILE,                                                         &
!     & RHOKM_U_1,RHOKM_V_1,                                            &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &
     & CDR10M_U,CDR10M_V,                                               &
     & GPP,NPP,RESP_P,GPP_FT,NPP_FT,                                    &
     & RESP_P_FT,RESP_S,RESP_S_TOT,RESP_W_FT,                           &
     ! Lestevens 23apr13
     & NPP_FT_ACC,RESP_W_FT_ACC,                                        &
     & WT_EXT_TILE,FSMC,                                                &
     & FLANDG_U,FLANDG_V,                                               &
     &            U_S_CAB,CH_CAB,CD_CAB,                                &
!sxy
     &           SNAGE_TILE,RTSOIL_TILE,                                &
     &           GFLUX_TILE,SGFLUX_TILE,                                &
     & F_ROOT, sin_theta_latitude,                                      &
! Lestevens Sept2012: CasaCNP variables 
     & CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
     & NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
      ! end step of experiment, this step, step width, processor num
     & endstep, timestep_number, mype                                   &
     & )

      Use bl_diags_mod, Only :                                          &
          strnewbldiag
      use cable_iday_mod
      IMPLICIT NONE

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
!EAK
!------------------------ nstypes.h ----------------------------------
!jhan:further renovation of ths file may be necessary params are dependent on dataset
!jhan: ALSO nstypes_cable.h should be unecessary nsoil/soil is only difference
      !--- Number of non-vegetation surface types
      Integer, Parameter :: NNVG  = 4

      !--- Number of plant functional types.
      Integer, Parameter :: NPFT  = 13
      
      !--- Number of surface types.
      Integer, Parameter :: NTYPE =17 
      
      !--- Index of the surface type 'Soil'
      !Integer, Parameter :: SOIL  = 16 
      !dhb599, 20110615: change made as per Peter Vohralik, item 1:
      Integer, Parameter :: SOIL  = 14

!--- Land surface types :
!--- original veg. tiles 
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!--- for testing these tiles are set = 1:5 
!     6 - Broadleaf Tree
!     7 - Needleleaf Tree
!     8 - C3 Grass
!     9 - C4 Grass
!    10 - Shrub
!--- for testing these tiles are set = 0
!    11 - 0 
!    11 - 0
!    11 - 0
!--- original non-veg tiles moved to these indices
!     14 - Urban
!     15 - Water
!     16 - Soil
!     17 - Ice



      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                             ! IN Number of X points?
     &,ROWS                                                             &
                             ! IN Number of Y points?
     &,OFF_X                                                            &
                             ! Size of small halo in i.
     &,OFF_Y                                                            &
                             ! Size of small halo in j.
     &,HALO_I                                                           &
                             ! Size of halo in i direction.
     &,HALO_J                                                           &
                             ! Size of halo in j direction.
     &,LAND_PTS                                                         &
                             ! IN No of land points being processed.
     &,NTILES                                                           &
                             ! IN Number of land tiles per land point.
     &,LAND_INDEX(LAND_PTS)                                             &
                             ! IN Index of land points.
     &,SM_LEVELS                                                        &
                                   
     &,TILE_INDEX(LAND_PTS,NTILES)                                      &
!                            ! IN Index of tile points.
     &,TILE_PTS(NTILES)                                                 &
                             ! IN Number of tile points.
     &,CAN_MODEL                                                        &
                             ! IN Switch for thermal vegetation canopy.
     &,FORMDRAG                                                         &
                             ! IN Switch for orographic form drag
     &,FD_stab_dep           ! IN Switch to implement stability
!                            !    dependence of orog form drag

      REAL                                                              &
     & BQ_1(ROW_LENGTH,ROWS)                                            &
                             ! IN A buoyancy parameter for lowest atm
!                            !    level ("beta-q twiddle").
     &,BT_1(ROW_LENGTH,ROWS)                                            &
                             ! IN A buoyancy parameter for lowest atm
!                            !    level ("beta-T twiddle").
     &,CANHC_TILE(LAND_PTS,NTILES)                                      &
!                            ! IN Areal heat capacity of canopy for
!                            !    land tiles (J/K/m2).
     &,CANOPY(LAND_PTS,NTILES)                                          &
!                            ! IN Surface water for land tiles
!                            !    (kg/m2).
     &,CATCH(LAND_PTS,NTILES)                                           &
                             ! IN Surface capacity (max. surface water)
!                            !    of land tiles (kg/m2).
     &,CATCH_SNOW(LAND_PTS)                                             &
                             ! IN Snow interception capacity of NLT
!                            !    tile (kg/m2).
     &,DZSOIL                                                           &
                             ! IN Soil or land-ice surface layer
!                            !    thickness (m).
     &,FLAKE(LAND_PTS,NTILES)                                           &
                             ! IN Lake fraction.
     &,GC(LAND_PTS,NTILES)                                              &
                             ! IN "Stomatal" conductance to evaporation
!                            !    for land tiles (m/s).
     &,HCON(LAND_PTS)                                                   &
                                   ! IN Soil thermal conductivity
!                                  !    (W/m/K).
     &,HCONS(LAND_PTS)                                                  &
                             ! IN Soil thermal conductivity including
!                            !    effects of water and ice (W/m/K).
     &,HO2R2_OROG(LAND_PTS)                                             &
                             ! IN Peak to trough height of unresolved
!                            !    orography divided by 2SQRT(2) (m).
     &,OROG_DRAG_PARAM                                                  &
!                            ! IN Drag coefficient for orographic
!                            !    form drag
     &,ICE_FRACT(ROW_LENGTH,ROWS)                                       &
!                            ! IN Fraction of gridbox which is sea-ice.
     &,FLAND(LAND_PTS)                                                  &
                             ! IN Land fraction on land tiles.
     &,FLANDG(ROW_LENGTH,ROWS)                                          &
!                            ! IN Land fraction on all tiles.

     &,SNOW_TILE(LAND_PTS,NTILES)                                       &
!                            ! IN Lying snow on land tiles (kg/m2).
     &,PSTAR(ROW_LENGTH,ROWS)                                           &
                             ! IN Surface pressure (Pascals).
     &,QW_1(ROW_LENGTH,ROWS)                                            &
                             ! IN Total water content of lowest
!                            !    atmospheric layer (kg per kg air).
     &,RADNET_SICE(ROW_LENGTH,ROWS)                                     &
!                            ! IN Sea-ice net surface radiation (W/m2)
     &,RADNET_TILE(LAND_PTS,NTILES)                                     &
!                            ! IN Land tile net surface radiation (W/m2)
     &,ANTHROP_HEAT(NTILES)                                             &
!                            ! IN Anthropogenic Urban heat source (W/m2)
     &,SIL_OROG_LAND(LAND_PTS)                                          &
                             ! IN Silhouette area of unresolved
!                            !    orography per unit horizontal area
     &,SMVCST(LAND_PTS)                                                 &
                             ! IN Volumetric saturation point
!                            !    - zero at land-ice points.
     &,TILE_FRAC(LAND_PTS,NTILES)                                       &
!                            ! IN Tile fractions.
     &,TIMESTEP                                                         &
                             ! IN Timestep in seconds for EPDT calc.
     &,TL_1(ROW_LENGTH,ROWS)                                            &
                             ! IN Liquid/frozen water temperature for
!                            !    lowest atmospheric layer (K).
     &,TI(ROW_LENGTH,ROWS)                                              &
                             ! IN Temperature of sea-ice surface layer
!                            !    (K)
     &,TS1(LAND_PTS)                                                    &
                             ! IN Temperature of top soil or land-ice
!                            !    layer (K)
     &,TSTAR_TILE(LAND_PTS,NTILES)                                      &
!                            ! IN Tile surface temperatures (K).
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
!                            ! IN Land mean surface temperature (K).
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
!                            ! IN Open sea surface temperature (K).
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
!                            ! IN Sea-ice surface temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS) !                                      &
!                            ! IN Mean sea surface temperature (K).
!dhb: seperate the continuation lines (too many...)
     REAL &
     & Z_LAND(ROW_LENGTH,ROWS)                                          &
!                            ! IN Land height (m).
     &, TSTAR(ROW_LENGTH,ROWS)                                          &
                              ! IN Gridbox Mean Surface Temperature (K)
     &,VFRAC_TILE(LAND_PTS,NTILES)                                      &
!                            ! IN Fractional canopy coverage for
!                            !    land tiles.
     &,VSHR_LAND(ROW_LENGTH,ROWS)                                       &
!                            ! IN Magnitude of land sfc-to-lowest-level
!                            !    wind shear
     &,VSHR_SSI(ROW_LENGTH,ROWS)                                        &
!                            ! IN Mag. of mean sea sfc-to-lowest-level
!                            !    wind shear
     &,ZH(ROW_LENGTH,ROWS)                                              &
                             ! IN Height above surface of top of
!                            !    boundary layer (metres).
     &,Z0_TILE(LAND_PTS,NTILES)                                         &
!                            ! IN Tile roughness lengths (m).
     &,Z1_UV(ROW_LENGTH,ROWS)                                           &
                             ! IN Height of lowest uv level (m).
     &,Z1_TQ(ROW_LENGTH,ROWS)                                           &
                             ! IN Height of lowest tq level (m).
!                            !    Note, if the grid used is staggered in
!                            !    the vertical, Z1_UV and Z1_TQ can be
!                            !    different.
     & ,Charnock   ! Charnock parameter for sea surface

      REAL                                                              &
     & Z0H_SCM(ROW_LENGTH,ROWS)                                         &
                               ! IN Namelist input z0h (if >0)
                               !    (if <=0 use Z0HSEA)
                               !    Used in SCM Configurations
     &,Z0M_SCM(ROW_LENGTH,ROWS)! IN Namelist input z0m (if >0)
                               !    (if <=0 use standard Z0MSEA)
                               !    Used in SCM Configurations

      LOGICAL                                                           &
     & LAND_MASK(ROW_LENGTH,ROWS)                                       &
!                            ! IN .TRUE. for land; .FALSE. elsewhere.
     &,SU10                                                             &
                             ! IN STASH flag for 10-metre W wind.
     &,SV10                                                             &
                             ! IN STASH flag for 10-metre S wind.
     &,SQ1P5                                                            &
                             ! IN STASH flag for 1.5-metre sp humidity.
     &,ST1P5                                                            &
                             ! IN STASH flag for 1.5-metre temperature.
     &,SFME                                                             &
                             ! IN STASH flag for wind mixing energy flux
     &,SZ0HEFF                                                          &
                             ! IN STASH flag for Z0H_EFF
     &,lq_mix_bl                                                        &
                             ! IN TRUE if mixing ratios used in
!                            !    boundary layer code
     &,LTIMER                                                           &
                             ! IN Logical for TIMER.
     &,L_DUST                                                           &
                              ! IN switch for mineral dust
     &,L_CTILE                                                          &
                              ! IN switch for coastal tiling
     &,L_CO2_INTERACTIVE                                                &
                              ! IN Switch for 3D CO2 field
     &,L_spec_z0              ! IN T if using prescribed
!                             !    sea surface roughness lengths
!
      REAL, Intent(IN) :: SeaSalinityFactor
!       Factor allowing for the effect of the salinity of
!       sea water on the evaporative flux.
      INTEGER, Intent(IN) :: ISeaZ0T
!       Switch for the treatment of the thermal roughness
!       length over the sea.
      INTEGER, Intent(IN) :: I_SCRN_T_DIAG
!       Method of diagnosing the screen temperature
!
      INTEGER                                                           &
     & COR_UST                                                          &
!                            ! IN Switch for friction velocity corr.n
     &,COR_MO_ITER           ! IN Switch for MO iteration correction

!  Modified (INOUT) variables.

      REAL                                                              &
     & Z0MSEA(ROW_LENGTH,ROWS)
!                            ! INOUT Sea-surface roughness length for
!                            !       momentum (m).  F617.

!  Output variables.
!
      REAL                                                              &
     & ALPHA1(LAND_PTS,NTILES)                                          &
!                            ! OUT Gradients of saturated specific
!                            !     humidity with respect to temperature
!                            !     between the bottom model layer and
!                            !     tile surface
     &,ALPHA1_SICE(ROW_LENGTH,ROWS)                                     &
!                            ! OUT ALPHA1 for sea-ice.
     &,ASHTF(ROW_LENGTH,ROWS)                                           &
                             ! OUT Coefficient to calculate surface heat
!                            !     flux into sea-ice (W/m2/K)
     &,ASHTF_TILE(LAND_PTS,NTILES)                                      &
!                            ! OUT Coefficient to calculate surface heat
!                            !     flux into land tiles (W/m2/K)
     &,CD(ROW_LENGTH,ROWS)                                              &
                             ! OUT Bulk transfer coefficient for
!                            !      momentum.
     &,CD_SSI(ROW_LENGTH,ROWS)                                          &
!                            ! OUT Bulk transfer coefficient for
!                            !      momentum over sea mean.
     &,CH(ROW_LENGTH,ROWS)                                              &
                             ! OUT Bulk transfer coefficient for heat
!                            !     and/or moisture.
     &,RECIP_L_MO_SEA(ROW_LENGTH,ROWS)                                  &
!                            ! OUT Reciprocal of the Monin-Obukhov
!                            !     length for sea/ice points (m^-1).
     &,CH_SSI(ROW_LENGTH,ROWS)                                          &
!                            ! OUT Bulk transfer coefficient for heat
!                            !    and/or moisture over sea mean.
     &,CDR10M(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)              &
!                            ! OUT Reqd for calculation of 10m wind
!                            !     (u & v).
!                            !     NBB: This is output on the UV-grid,
!                            !     but with the first and last rows set
!                            !     to a "missing data indicator".
!                            !     Sea-ice leads ignored.
     &,CHR1P5M(LAND_PTS,NTILES)                                         &
!                            ! OUT Reqd for calculation of 1.5m temp for
!                            !     land tiles.
     &,CHR1P5M_SICE(ROW_LENGTH,ROWS)                                    &
!                            ! OUT CHR1P5M for sea and sea-ice
!                            !     (leads ignored).
     &,E_SEA(ROW_LENGTH,ROWS)                                           &
                             ! OUT Evaporation from sea times leads
!                            !     fraction (kg/m2/s). Zero over land.
     &,FME(ROW_LENGTH,ROWS)                                             &
                             ! OUT Wind mixing energy flux (Watts/sq m).
     &,FQW_1(ROW_LENGTH,ROWS)                                           &
                             ! OUT "Explicit" surface flux of QW (i.e.
!                            !     evaporation), on P-grid (kg/m2/s).
!                            !     for whole grid-box
     &,FQW_TILE(LAND_PTS,NTILES)                                        &
!                            ! OUT Local FQW_1 for land tiles.
     &,FQW_ICE(ROW_LENGTH,ROWS)                                         &
!                            ! OUT GBM FQW_1 for sea-ice.
     &,FTL_1(ROW_LENGTH,ROWS)                                           &
                             ! OUT "Explicit" surface flux of TL = H/CP.
!                            !     (sensible heat / CP). grid-box mean
     &,FTL_TILE(LAND_PTS,NTILES)                                        &
!                            ! OUT Local FTL_1 for land tiles.
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
!                            ! OUT GBM FTL_1 for sea-ice.
     &,FRACA(LAND_PTS,NTILES)                                           &
                             ! OUT Fraction of surface moisture flux
!                            !     with only aerodynamic resistance
!                            !     for land tiles.
     &,H_BLEND_OROG(ROW_LENGTH,ROWS)                                    &
!                            ! OUT Blending height for orographic
!                            !     roughness
     &,H_SEA(ROW_LENGTH,ROWS)                                           &
                             ! OUT Surface sensible heat flux over sea
!                            !     times leads fraction (W/m2).
!                            !     Zero over land.
     &,RHOSTAR(ROW_LENGTH,ROWS)                                         &
!                            ! OUT Surface air density
     &,RESFS(LAND_PTS,NTILES)                                           &
!                            ! OUT Combined soil, stomatal and
!                            !     aerodynamic resistance factor for
!                            !     fraction 1-FRACA of land tiles
     &,RESFT(LAND_PTS,NTILES)                                           &
!                            ! OUT Total resistance factor
!                            !     FRACA+(1-FRACA)*RESFS for snow-free
!                            !     tiles, 1 for snow and land-ice.
     &,RIB(ROW_LENGTH,ROWS)                                             &
                             ! OUT Mean bulk Richardson number for
!                            !     lowest layer
     &,RIB_TILE(LAND_PTS,NTILES)                                        &
!                            ! OUT RIB for land tiles.
     &,FB_SURF(ROW_LENGTH,ROWS)                                         &
!                            ! OUT Surface flux buoyancy over
!                            !     density (m^2/s^3)
     &,U_S(ROW_LENGTH,ROWS)                                             &
                             ! OUT Surface friction velocity (m/s)
     &,Q1_SD(ROW_LENGTH,ROWS)                                           &
                             ! OUT Standard deviation of turbulent
!                            !     fluctuations of surface layer
!                            !     specific humidity (kg/kg).
     &,T1_SD(ROW_LENGTH,ROWS) !                                          &
                             ! OUT Standard deviation of turbulent
!                            !     fluctuations of surface layer
!                            !     temperature (K).
!dhb: seperate continuation lines (too many)...
     real &
     & Z0HSSI(ROW_LENGTH,ROWS)                                          &
                             ! OUT Roughness length for heat and
!                            !     moisture over sea/sea ice (m)
     &,Z0H(LAND_PTS)                                        &
     &,Z0H_TILE(LAND_PTS,NTILES)                                        &
!                            ! OUT Tile roughness lengths for heat
!                            !     and moisture
     &,Z0H_EFF(ROW_LENGTH,ROWS)                                         &
!                            ! OUT Effective roughness length for
!                            !     heat, moisture (m)
     &,Z0M_GB(ROW_LENGTH,ROWS)                                          &
                             ! OUT Gridbox mean roughness length 
!                            !     for momentum (m)
     &,Z0MSSI(ROW_LENGTH,ROWS)                                          &
                             ! OUT Roughness length for momentum
!                            !     over sea/sea ice (m)
     &,Z0M(LAND_PTS)                                        &
     &,Z0M_TILE(LAND_PTS,NTILES)                                        &
!                            ! OUT Tile roughness lengths for momentum
     &,Z0M_EFF(ROW_LENGTH,ROWS)                                         &
!                            ! OUT Effective roughness length for
!                            !     momentum
     &,RHO_ARESIST(ROW_LENGTH,ROWS)                                     &
!                            ! OUT RHOSTAR*CD_STD*VSHR  for SCYCLE
     &,ARESIST(ROW_LENGTH,ROWS)                                         &
!                            ! OUT 1/(CD_STD*VSHR)      for SCYCLE
     &,RESIST_B(ROW_LENGTH,ROWS)                                        &
!                            ! OUT (1/CH-1/CD_STD)/VSHR for SCYCLE
     &,RHO_ARESIST_TILE(LAND_PTS,NTILES)                                &
!                            ! OUT RHOSTAR*CD_STD*VSHR on land tiles
     &,ARESIST_TILE(LAND_PTS,NTILES)                                    &
!                            ! OUT 1/(CD_STD*VSHR) on land tiles
     &,RESIST_B_TILE(LAND_PTS,NTILES)                                   &
!                            ! OUT (1/CH-1/CD_STD)/VSHR on land tiles
     &,R_B_DUST(ROW_LENGTH,ROWS,NDIV)                                   &
!                            !OUT surf layer res for dust
     &,CD_STD_DUST(ROW_LENGTH,ROWS)                                     &
!                            !OUT Bulk transfer coef. for
!                            ! momentum, excluding orographic effects
     &,U_S_STD_TILE(LAND_PTS,NTILES)
!                            ! OUT Surface layer scaling velocity
!                            ! for tiles excluding orographic
!                            ! form drag (m/s).
      Type (Strnewbldiag) :: BL_diag


! Surface exchange coefficients;passed to subroutine IMPL_CAL
      REAL                                                              &
     & RHO_CD_MODV1(ROW_LENGTH,ROWS)                                    &
!                            ! OUT rhostar*cD*vshr before horizontal
!                            !     interpolation output as a diagnostic.
     &,RHOKH_1(LAND_PTS,NTILES)                                         &
!                            ! OUT Surface exchange coefficient for land
!                            !     tiles.
     &,RHOKH_1_SICE(ROW_LENGTH,ROWS)                                    &
!                            ! OUT Surface exchange coefficient for sea
!                            !     or sea-ice.
     &,RHOKM_1(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)             &
!                            ! OUT For momentum. NB: This is output on
!                            !     UV-grid, but with the first and last
!                            !     rows set to "missing data indicator".
     &,RHOKM_LAND(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)          &
!                            ! OUT For land momentum. NB: This is output
!                            !     on UV-grid, but with the first and
!                            !      last rows set to "missing data".
     &,RHOKM_SSI(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)           &
!                            ! OUT For mean sea mom. NB: This is output
!                            !     on UV-grid, but with the first and
!                            !     last rows set to "missing data".
     &,RHOKPM(LAND_PTS,NTILES)                                          &
!                            ! OUT Mixing coefficient for land tiles.
     &,EPOT_TILE(LAND_PTS,NTILES)                                       &
                                   ! OUT EPOT for land tiles.
     &,RHOKPM_POT(LAND_PTS,NTILES)                                      &
!                                  ! OUT Potential evaporation
!                                  !     exchange coeff.
     &,RHOKPM_SICE(ROW_LENGTH,ROWS)
!                            ! OUT Mixing coefficient for sea-ice.

!------------------------------------------------------
!     EAK

      INTEGER                                                           &
     &  N_ROWS                                                           &
!                  ! Local number of rows in a v field
                                   ! IN No. of soil moisture levels
     &,CO2_DIM_LEN                                                      &
                                   ! IN Length of a CO2 field row.
     &,CO2_DIM_ROW                                                      &

     &, day                                                             &
     &, total_nsteps                                                    &
                                ! Total number of steps in run
     &, SOIL_TYPE(row_length,rows)                                      &
     &, VEG_TYPE(row_length,rows)                                       &

     &, ISNOW_FLG3L(LAND_PTS,NTILES)                             &
     &, DIM_CS1, DIM_CS2                                                 &
                                ! soil carbon dimensions
     &, LAND_PTS_TRIF                                                    &
!                 ! IN For dimensioning land fields
     &, NPFT_TRIF                                                        
                  ! IN For dimensioning PFT fields available only
!                 !    with TRIFFID. Set to NPFT when TRIFFID on,
!                 !    set to 1 when TRIFFID off.
      Logical                                                           &
     & l_cable
      LOGICAL, DIMENSION(LAND_PTS,NTILES) :: L_TILE_PTS

! 
      REAL                                                       &
     &  alb_tile(land_pts,ntiles,4)                              &
     &, albsoil(land_pts)                                        &
     &,CO2_MMR                                                          &
                                   ! IN CO2 Mass Mixing Ratio
     &,CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)                                  &
!                                  ! IN 3D CO2 field if required.
     &, SNOW_DEPTH3L(LAND_PTS,NTILES,3)                          &
     &, SNOW_MASS3L(LAND_PTS,NTILES,3)                           &
     &, SNOW_COND(LAND_PTS,NTILES,3)                             &
     &, SNOW_TMP3L(LAND_PTS,NTILES,3)                            &
     &, SNOW_RHO3L(LAND_PTS,NTILES,3)                            &
     &, SNOW_RHO1L(LAND_PTS,NTILES)                              &
     &, SNAGE_TILE(LAND_PTS,NTILES)                              &
     &, SMCL_TILE(LAND_PTS,NTILES,sm_levels)                     &
     &, STHU_TILE(LAND_PTS,NTILES,sm_levels)                     &
     &, STHF_TILE(LAND_PTS,NTILES,sm_levels)                     &
     &, TSOIL_TILE(LAND_PTS,NTILES,sm_levels)                    &
     &, T_SURF_TILE(LAND_PTS,NTILES)                             &
     &, RTSOIL_TILE(LAND_PTS,NTILES)                             &
     &, GFLUX_TILE(LAND_PTS,NTILES)                              &
     &, SGFLUX_TILE(LAND_PTS,NTILES)                             &
     &, BEXP(LAND_PTS)                                           &
     &, SATCON(LAND_PTS)                                         &
     &, SATHH(LAND_PTS)                                          &
     &, SMVCCL(LAND_PTS)                                         &
     &, SMVCWT(LAND_PTS)                                         & 
     &, HCAP(LAND_PTS)                                           &
! Lestevens Sept2012: CasaCNP variables 
     &, CPOOL_TILE(LAND_PTS,NTILES,10)                           &
     &, NPOOL_TILE(LAND_PTS,NTILES,10)                           &
     &, PPOOL_TILE(LAND_PTS,NTILES,12)                           &
     &, SOIL_ORDER(LAND_PTS)                                     &
     &, NIDEP(land_pts)                                          &
     &, NIFIX(land_pts)                                          &
     &, PWEA(land_pts)                                           &
     &, PDUST(land_pts)                                          &
     &, GLAI(land_pts,NTILES)                                    &
     &, PHENPHASE(land_pts,NTILES)                               &

     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)                                     &
     &, surf_down_sw(row_length,rows,4)                                 &
     &, cos_zenith_angle(row_length,rows)                               &
     &, sin_theta_latitude(row_length,rows)                             &
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, LAT(ROW_LENGTH,ROWS)                                            &
     &, LONG(ROW_LENGTH,ROWS)                                           &
     &, time_sec                                                        &
     &, DS_RATIO                                                        &
     &, LH                                                              &
                           ! Latent heat (J/K/kg).
     &, D_T                                                             &

     &,FTL_TILE_CAB(LAND_PTS,NTILES)                                 &
     &,FTL_CAB(LAND_PTS)                                             &
     &,LE_TILE_CAB(LAND_PTS,NTILES)                                  &
     &,LE_CAB(LAND_PTS)                                              &
     &,TSTAR_TILE_CAB(LAND_PTS,NTILES)                               &
     &,TSTAR_CAB(LAND_PTS)                                           &
     &,SMCL_CAB(LAND_PTS,SM_LEVELS)                          &
     &,TSOIL_CAB(LAND_PTS,SM_LEVELS)                         &
     &,USTAR_CAB(LAND_PTS)                                           &
!sxy
     &,CH_CAB(LAND_PTS)                                              &
     &,CD_CAB(LAND_PTS)                                              &
     &,U_S_CAB(LAND_PTS)                                             &
     &,SURF_HTF_CAB(LAND_PTS)                                         


!-------------------------------------------------------------------------
      REAL                                                              &
     & CDR10M_U(ROW_LENGTH,ROWS)                                        &
                                   ! OUT Ratio of CD's reqd for
!                                  !     calculation of 10 m wind. On
!                                  !     U-grid; comments as per RHOKM.
     &,CDR10M_V(ROW_LENGTH,ROWS)                                      &
                                   ! OUT Ratio of CD's reqd for
!                                  !     calculation of 10 m wind. On
!                                  !     V-grid; comments as per RHOKM.
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
                                   ! OUT Land frac (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
!                                  !     land tiles (W/m2)
     &,FLANDG_V(ROW_LENGTH,ROWS)                                      &
                                   ! OUT Land frac (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,FSMC(LAND_PTS,NPFT)                                              &
                                       ! OUT Moisture availability factor.
     &,GS(LAND_PTS)                                                     &
                                   ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
     &,LW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface downward LW radiation
!                                  !    (W/m2).
     &,SW_DOWN(ROW_LENGTH,ROWS)                                         &  
                                   ! Surface downward SW radiation (W/m2).

     &,SW_TILE(LAND_PTS,NTILES)                                         &
                                   ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
     &,LE_TILE(LAND_PTS,NTILES)                                         &
                                   ! OUT Surface latent heat flux for
!                                  !     land tiles
     &,U(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)                 &
!                                  ! IN W'ly wind component (m/s)
     &,V(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:N_ROWS+OFF_Y)               &
!                                  ! IN S'ly wind component (m/s)
     &,U_P(ROW_LENGTH,ROWS)                                           &
                                   ! IN U_1 on P-grid.
     &,V_P(ROW_LENGTH,ROWS)                                           &
                                   ! IN V_1 on P-grid.

     &,STHF(LAND_PTS,SM_LEVELS)                                         &
                                   ! IN Frozen soil moisture content of
!                                  !    each layer as a fraction of
!                                  !    saturation.
     &,STHU(LAND_PTS,SM_LEVELS)                                         &
                                   ! IN Unfrozen soil moisture content
!                                  !    of each layer as a fraction of
!                                  !    saturation.
     &,SOIL_LAYER_MOISTURE(LAND_PTS,SM_LEVELS)                         &
                                               !IN soil moisture
!                 ! per layer (kg m-2)
     &,T_SOIL(LAND_PTS,SM_LEVELS)                                       &
                                   ! IN Soil temperatures (K).
     &,CS(LAND_PTS,DIM_CS1)                                             &
                              ! IN Soil carbon (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                   ! IN Fractions of surface types.
     &,CANHT_FT(LAND_PTS,NPFT)                                          &
                                   ! IN Canopy height (m)
     &,PHOTOSYNTH_ACT_RAD(ROW_LENGTH,ROWS)                              &
!                                  ! IN Net downward shortwave radiation
!                                  !    in band 1 (w/m2).
     &,LAI_FT(LAND_PTS,NPFT)                                            &
                                   ! IN Leaf area index
     !kdcorbin, 11/10 - changed from NPFT
     &,G_LEAF(LAND_PTS,NTILES)                                          &
                                   ! OUT Leaf turnover rate (/360days).
     &,G_LEAF_ACC(LAND_PTS,NPFT)                                        &
                                   ! INOUT Accumulated G_LEAF
! Lestevens 17may13 - change npft to ntiles
     &,NPP_FT_ACC(LAND_PTS_TRIF,NTILES)                                 &
!    &,NPP_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                              &
!                                  ! INOUT Accumulated NPP_FT
     &,RESP_W_FT_ACC(LAND_PTS_TRIF,NTILES)                              &
!    &,RESP_W_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                           &
!                                  ! INOUT Accum RESP_W_FT
     &,RESP_S_ACC(LAND_PTS_TRIF,DIM_CS1)                                &
                                   ! INOUT Accumulated RESP_S     

!     &,RHOKM_U_1(ROW_LENGTH,ROWS)                                       &
!                                   ! OUT Exchange coefficients for
!!                                  !     momentum (on U-grid, with 1st
!!                                  !     and last rows undefined or, at
!!                                  !     present, set to "missing data")
!     &,RHOKM_V_1(ROW_LENGTH,N_ROWS)                                     &
!                                   ! OUT Exchange coefficients for
!!                                  !     momentum (on V-grid, with 1st
!!                                  !     and last rows undefined or, at
!!                                  !     present, set to "missing data")
!     &,TAUX_1(ROW_LENGTH,ROWS)                                          &
!                                   ! OUT W'ly component of surface wind
!     &,TAUY_1(ROW_LENGTH,N_ROWS)                                        &
!                                   ! OUT S'ly component of surface wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       &
                                   ! OUT W'ly component of sfc wind
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
     &,TAUX_SSI(ROW_LENGTH,ROWS)                                        &
                                   ! OUT W'ly component of sfc wind
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     &
                                   ! OUT S'ly component of sfc wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
     &,TAUY_SSI(ROW_LENGTH,N_ROWS)                                      &
                                   ! OUT S'ly component of sfc wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
!                                  !     set to missing data
     &,GPP(LAND_PTS)                                                    &
                                   ! OUT Gross primary productivity
!                                  !     (kg C/m2/s).
     &,NPP(LAND_PTS)                                                    &
                                   ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
     &,RESP_P(LAND_PTS)                                                 &
                                   ! OUT Plant respiration (kg C/m2/s).
!     &,G_LEAF(LAND_PTS,NPFT)                                            &
!                                   ! OUT Leaf turnover rate (/360days).
     !kdcorbin, 11/10 - changed from NPFT
     &,GPP_FT(LAND_PTS,NTILES)                                          &
                                   ! OUT Gross primary productivity
!                                  !     on PFTs (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,NPP_FT(LAND_PTS,NTILES)                                          &
                                   ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,RESP_P_FT(LAND_PTS,NTILES)                                       &
                                   ! OUT Plant respiration on PFTs
!                                  !     (kg C/m2/s).
     &,RESP_S(LAND_PTS,DIM_CS1)                                         &
                                ! OUT Soil respiration (kg C/m2/s).
     &,RESP_S_TOT(DIM_CS2)                                              &
                                  ! OUT Total soil respiration
                                  ! (kg C/m2/s).
     &,RESP_W_FT(LAND_PTS,NPFT)                                         &
                                   ! OUT Wood maintenance respiration
!                                  !     (kg C/m2/s).
     &,WT_EXT_TILE(LAND_PTS,SM_LEVELS,NTILES)                           &
!                                  ! IN Fraction of evapotranspiration
!                                  !    which is extracted from each
!                                  !    soil layer by each tile.
     &, F_ROOT(sm_levels)

!  Symbolic constants ------------------------------------------------

!   (1) UM-wide common parameters.

!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
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
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
! C_PERMA start

      ! Specific heat capacity of water vapour (J/kg/K)
      REAL,PARAMETER:: HCAPV=1850.0

      ! Specific heat capacity of water (J/kg/K)
      REAL,PARAMETER:: HCAPW=4180.0

      ! Specific heat capacity of ice (J/kg/K)
      REAL,PARAMETER:: HCAPI=2100.0

      ! Density of ice (kg/m3)
      REAL,PARAMETER:: RHO_ICE=917

      ! Rate of change of ice potential with temperature
      ! RHO_ICE*LF/ZERODEGC*1/(RHO_WATER*G) (m/K)
      REAL,PARAMETER:: DPSIDT=114.3

! C_PERMA end
!EAK
! C_SOILH start
      ! No. of soil layers (must = NSOIL).
      REAL,PARAMETER:: PSOIL=4

      ! Tunable characteristic freq (rad/s)
      REAL,PARAMETER:: OMEGA1=3.55088E-4

      ! Density of lying snow (kg per m**3)
      REAL,PARAMETER:: RHO_SNOW=250.0

      ! Depth of `effective' snow surface layer (m)
      REAL,PARAMETER:: DEFF_SNOW=0.1

      ! Thermal conductivity of lying snow (Watts per m per K).
      REAL,PARAMETER:: SNOW_HCON=0.265

      ! Thermal capacity of lying snow (J/K/m3)
      REAL,PARAMETER:: SNOW_HCAP=0.63E6

! C_SOILH end
! CSIGMA start
      ! Stefan-Boltzmann constant (W/m**2/K**4).
      REAL, PARAMETER ::  SBCON=5.67E-8
! CSIGMA end

! Derived local parameters.
      REAL GRCP,LS
      PARAMETER (                                                       &
     & GRCP=G/CP                                                        &
     &,LS=LF+LC                                                         &
                           ! Latent heat of sublimation.
     & )

!   (2) Boundary Layer local parameters.

! Description:
!   This deck sets up the parameter LB
!
! Current Code Owner: Z Gardner
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 5.3      16/08/01 Added in header and changed definitions. Z Gardner
! 5.5      17/04/03 Remove reference to obsolete section
!                   A03_7A. T.White
!
! Declarations:
! Start blend_h
! Description:
!   This file sets the value of the variable LB
!
! Current Code Owner:
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.3   25/09/01  Portability changes.  Z. Gardner

      REAL,PARAMETER:: LB = 20.0 ! Blending height (m).
! End Blend_h
! C_DENSTY for subroutine SF_EXCH
      REAL,PARAMETER:: RHOSEA = 1026.0 ! density of sea water (kg/m3)
      REAL,PARAMETER:: RHO_WATER = 1000.0! density of pure water (kg/m3)
! C_DENSTY end
!!----------------------------------------------------------------------
!!!-----------COMDECK C_ROUGH FOR SUBROUTINE SF_EXCH----------
! Sea ice parameters
! Z0FSEA = roughness length for free convective heat and moisture
!          transport over the sea (m).
!          DUMMY VARIABLE - Only used in 7A boundary layer scheme
! Z0HSEA = roughness length for heat and moisture transport
!          over the sea (m).
! Z0MIZ  = roughness length for heat, moisture and momentum over
!          the Marginal Ice Zone (m).
! Z0SICE = roughness length for heat, moisture and momentum over
!          sea-ice (m).
      REAL :: Z0HSEA
      REAL :: Z0MIZ
      REAL :: Z0SICE

      COMMON  /RUN_BLICE/Z0HSEA,Z0MIZ,Z0SICE

!!----------------------------------------------------------------------
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
      REAL                                                              &
     & Z0H_Z0M(17)                 ! Ratio of roughness length for heat
!                                ! to roughness length for momentum.
!----------------------------------------------------------------------
!                        BT   NT   C3G  C4G  Shr  Urb  Wat  Soil Ice
!----------------------------------------------------------------------
      DATA Z0H_Z0M  /.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1/

      REAL H_BLEND_MIN
      PARAMETER (                                                       &
     & H_BLEND_MIN=0.0                                                  &
                             ! Minimum blending height.
     &)

! Start blopt8a

! Description:
!   Permissible settings for BL options.
!
! Current Code Owner: J. M. Edwards
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 6.2      27/01/06 Original code.  J. M. Edwards
!
      INTEGER, PARAMETER :: Off = 0  ! Switch disabled
      INTEGER, PARAMETER :: On  = 1  ! Switch enabled
!
      INTEGER, PARAMETER :: BrownGrant97 = 1
      INTEGER, PARAMETER :: BrownGrant97_limited = 2
!     Options for non-gradient stress following
!     Brown and Grant (1997), version 2 including a limit on its size
!
!     Options for flux gradient formulation
      INTEGER, PARAMETER :: Locketal2000   = 0
!       Flux gradients as in Lock et al. (2000)
      INTEGER, PARAMETER :: HoltBov1993 = 1
!       Flux gradients as in Lock et al (2000) but using
!       coefficients from Holtslag and Boville (1993)
      INTEGER, PARAMETER :: LockWhelan2006 = 2
!       Flux gradients as in Lock and Whelan (2006)
!
!     Options for form drag
      INTEGER, PARAMETER :: No_drag         = 0
      INTEGER, PARAMETER :: Effective_z0    = 1
      INTEGER, PARAMETER :: Explicit_stress = 2
!
!     Options for marine boundary layers
      INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
      INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory
      INTEGER, PARAMETER :: DynDiag_ZL = 1
!       The ratio of the height of the inversion to the surface
!       Obukhov length is used a dynamic criterion in the
!       diagnosis of BL types
!
      INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
!
!     Options for stable boundary layers
      INTEGER, PARAMETER ::  Long_tails           = 0
      INTEGER, PARAMETER ::  Sharpest             = 1
      INTEGER, PARAMETER ::  Sharp_sea_long_land  = 2
      INTEGER, PARAMETER ::  Mes_tails            = 3
      INTEGER, PARAMETER ::  Louis_tails          = 4
      INTEGER, PARAMETER ::  Depth_based          = 5
      INTEGER, PARAMETER ::  Sharp_sea_mes_land   = 6
      INTEGER, PARAMETER ::  LEM_stability        = 7
      INTEGER, PARAMETER ::  Sharp_sea_Louis_land = 8

!     Options for Prandtl number (in local Ri scheme)
      INTEGER, PARAMETER ::  Constant_SBL = 0
      INTEGER, PARAMETER ::  LockMailhot2004 = 1

! End blopt8a

!   External subprograms called.

      EXTERNAL QSAT_mix,SF_OROG,SF_RESIST,SF_RIB_SEA,                   &
     & SF_RIB_LAND,FCDCH_SEA,FCDCH_LAND,SF_FLUX_SEA,                    &
     & SF_FLUX_LAND,STDEV1_SEA,STDEV1_LAND,SF_OROG_GB,                  &
     & SFL_INT_SEA,SFL_INT_LAND                                         &
     & ,DUSTRESB
      EXTERNAL TIMER


!   Define local storage.

!   (a) Workspace.

      REAL                                                              &
     & QS1(ROW_LENGTH,ROWS)        ! Sat. specific humidity
!                                  ! qsat(TL_1,PSTAR)

!  Workspace for sea and sea-ice leads
      REAL                                                              &
     & CD_SEA(ROW_LENGTH,ROWS)                                          &
                                   ! Drag coefficient
     &,CH_SEA(ROW_LENGTH,ROWS)                                          &
                                   ! Transfer coefficient for heat and
!                                  ! moisture
     &,QSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! Surface saturated sp humidity
     &,RIB_SEA(ROW_LENGTH,ROWS)                                         &
                                   ! Bulk Richardson number
     &,Z0H_SEA(ROW_LENGTH,ROWS)                                         &
                                   ! Roughness length for heat and
!                                  ! moisture transport
     &,Z0M_SEA(ROW_LENGTH,ROWS)                                         &
                                   ! Open sea roughness length for
!                                  ! momentum transport.
     &,DB_SEA(ROW_LENGTH,ROWS)                                          &
                                   ! Buoyancy difference for sea points
     &,V_S_SEA(ROW_LENGTH,ROWS)    ! Surface layer scaling velocity
!                                  ! for sea points (m/s).

!  Workspace for sea-ice and marginal ice zone
      REAL                                                              &
     & CD_ICE(ROW_LENGTH,ROWS)                                          &
                                   ! Drag coefficient
     &,CD_LAND(ROW_LENGTH,ROWS)                                         &
                                   ! Bulk transfer coefficient for
!                                  !      momentum over land.
     &,CD_MIZ(ROW_LENGTH,ROWS)                                          &
                                   ! Drag coefficient
     &,CH_ICE(ROW_LENGTH,ROWS)                                          &
                                   ! Transfer coefficient for heat and
!                                  ! moisture
     &,CH_MIZ(ROW_LENGTH,ROWS)                                          &
                                   ! Transfer coefficient for heat and
!                                  ! moisture
     &,QSTAR_ICE(ROW_LENGTH,ROWS)                                       &
                                   ! Surface saturated sp humidity
     &,RIB_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! Bulk Richardson number
     &,RIB_MIZ(ROW_LENGTH,ROWS)                                         &
                                   ! Bulk Richardson number
     &,Z0_ICE(ROW_LENGTH,ROWS)                                          &
                                   ! Roughness length.
     &,Z0_MIZ(ROW_LENGTH,ROWS)                                          &
                                   ! Roughness length.
     &,DB_ICE(ROW_LENGTH,ROWS)                                          &
                                   ! Buoyancy difference for sea ice
     &,V_S_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! Surface layer scaling velocity
!                                  ! for sea ice (m/s).
     &,V_S_MIZ(ROW_LENGTH,ROWS)                                         &
                                   ! Surface layer scaling velocity
!                                  ! for marginal sea ice (m/s).
     &,RECIP_L_MO_ICE(ROW_LENGTH,ROWS)                                  &
!                                  ! Reciprocal of the Monin-Obukhov
!                                  ! length for sea ice (m^-1).
     &,RECIP_L_MO_MIZ(ROW_LENGTH,ROWS)                                  &
!                                  ! Reciprocal of the Monin-Obukhov
!                                  ! length for marginal sea ice (m^-1).
     &,RHO_ARESIST_LAND(ROW_LENGTH,ROWS)
!                            ! Land mean of rho_aresist_tile
      REAL                                                              &
     & Z1_TQ_SEA(ROW_LENGTH,ROWS)
!                            ! Height of lowest model level
!                            ! relative to sea.
      INTEGER                                                           &
     & SICE_INDEX(ROW_LENGTH*ROWS,2)                                    &
!                                  ! Index of sea-ice points
     &,NSICE                       ! Number of sea-ice points.

!  Workspace for land tiles
      REAL                                                              &
     & CD_STD(LAND_PTS,NTILES)                                          &
                                   ! Local drag coefficient for calc
!                                  ! of interpolation coefficient
     &,CD_TILE(LAND_PTS,NTILES)                                         &
                                   ! Drag coefficient
     &,CH_TILE(LAND_PTS,NTILES)                                         &
                                   ! Transfer coefficient for heat and
!                                  ! moisture
     &,CHN(LAND_PTS)                                                    &
                                   ! Neutral value of CH.
     &,CE_D_SH(LAND_PTS)                                                &
                                   ! Exposure coefficient * diffusivity
!                                  ! of water vapour * Sherwood number
!                                  ! for canopy snow (m^2/s).
     &,DQ(LAND_PTS)                                                     &
                                   ! Sp humidity difference between
!                                  ! surface and lowest atmospheric lev
     &,EPDT(LAND_PTS)                                                   &
                                   ! "Potential" Evaporation * Timestep
     &,FZ0(LAND_PTS)                                                    &
                                   ! Aggregation function for Z0.
     &,PSTAR_LAND(LAND_PTS)                                             &
                                   ! Surface pressure for land points.
     &,QSTAR_TILE(LAND_PTS,NTILES)                                      &
                                   !Surface saturated sp humidity.
     &,RH_CAN(LAND_PTS)                                                 &
                                   ! Canopy air to surface resistance
!                                  ! (s/m).
     &,RHOKM_1_TILE(LAND_PTS,NTILES)                                    &
!                                  ! Momentum exchange coefficient.
     &,WIND_PROFILE_FACTOR(LAND_PTS,NTILES)                             &
!                                  ! For transforming effective surface
!                                  ! transfer coefficients to those
!                                  ! excluding form drag.
     &,Z0M_GB_LAND(LAND_PTS)                                            &
                                   ! GBM momentum land roughness length 
     &,Z0H_GB_LAND(LAND_PTS)                                            &
                                   ! GBM land roughness length for heat
     &,Z0M_EFF_TILE(LAND_PTS,NTILES)                                    &
!                                  ! Effective momentum roughness length
     &,DB_TILE(LAND_PTS,NTILES)                                         &
                                   ! Buoyancy difference for surface
!                                  ! tile
     &,V_S_TILE(LAND_PTS,NTILES)                                        &
                                   ! Surface layer scaling velocity
!                                  ! for tiles (m/s).
     &,V_S_STD(LAND_PTS,NTILES)                                         &
!                                  ! Surface layer scaling velocity
!                                  ! for tiles excluding orographic
!                                  ! form drag (m/s).
     &,U_S_ITER_TILE(LAND_PTS,NTILES)                                   &
!                                  ! Scaling velocity from middle of 
!                                  ! MO scheme - picked up in error by 
!                                  ! dust code!
     &,VSHR(ROW_LENGTH,ROWS)                                            &
!                                  ! Level 1 to surface wind shear
     &,RECIP_L_MO_TILE(LAND_PTS,NTILES)
!                                  ! Reciprocal of the Monin-Obukhov
!                                  ! length for tiles (m^-1).

!   (b) Scalars.

      INTEGER                                                           &
     & I,J                                                              &
                   ! Loop counter (horizontal field index).
     &,K                                                                &
                   ! Loop counter (tile field index).
     &,L                                                                &
                   ! Loop counter (land point field index).
     &,N                                                                &
                   ! Loop counter (tile index).
     &,IDIV                                                             &
                   ! Loop counter (dust division).
     &,JITS        ! Counter for iteration for Z0H
      REAL                                                              &
     & TAU                                                              &
                   ! Magnitude of surface wind stress over sea.
     &,ZETAM                                                            &
                   ! Temporary in calculation of CHN.
     &,ZETAH                                                            &
                   ! Temporary in calculation of CHN.
     &,ZETA1                                                            &
                   ! Work space
     &,Z0                                                               &
                   ! yet more workspace
     &,USTR_L                                                           &
                   ! Low-wind estimate of friction velocity
     &,USTR_N                                                           &
                   ! Neutral estimate of friction velocity
     &,TOL_USTR_N                                                       &
                   ! Tolerance for USTR_N
     &,TOL_USTR_L                                                       &
                   ! Tolerance for USTR_L (see below)
     &,TINY_SNOW   ! Smallest amount of snow for > test, to avoid NaNs
!EAK
       INTEGER, PARAMETER  :: ntest = 0 !  for prints

      ! end step of experiment, step width, processor num
      integer :: endstep, timestep_number, mype

!	YMA: newly defined variables for the changed z0m and z0h expressions
!	1). variables in the Charnock expression in z0m form 
!        REAL CK4Z0M, CK4Z0ML, CK4Z0MH, XCK4Z0M
!        REAL UMD4CK, RT4CK
!dhb:
	real :: CK4Z0M=0.
	real :: CK4Z0ML=0.
	real :: CK4Z0MH=0.
	real :: XCK4Z0M=0.
	real :: UMD4CK=0.
	real :: RT4CK=0.
!	2). variables in the z0h form
!	REAL RT4Z0H, ALGRRDIF4Z0H, SH4Z0H, ALGZ0T4Z0H
!dhb:
	real :: RT4Z0H=0.
	real :: ALGRRDIF4Z0H=0.
	real :: SH4Z0H=0.
	real :: ALGZ0T4Z0H=0.
!       3.)  define variables for temp uses
	REAL :: EXPOSI=0.
	REAL :: EXNEGA=0.

!	YMA: define variables in the Charnock expression
        CK4Z0ML = 0.011
        CK4Z0MH = 0.032
        UMD4CK = 15.
	RT4CK = 0.2

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFEXCH  ',3)
      ENDIF

!-----------------------------------------------------------------------
!!  0. Initialise FTL_TILE and RIB_TILE on all tiles at all points,
!!     to allow STASH to process these as diagnostics.
!-----------------------------------------------------------------------
      DO N=1,NTILES
!CDIR NODEP
        DO L=1,LAND_PTS
          FTL_TILE(L,N) = 0.0
          RIB_TILE(L,N) = 0.0
          Z0M_TILE(L,N) = 0.0
          U_S_STD_TILE(L,N)=0.0
!EAK
          IF(l_cable .and. TILE_FRAC(L,N) .eq. 0.0 ) THEN
            VFRAC_TILE(L,N) = 0.0
            CANHC_TILE(L,N) = 0.0
            CD_TILE(L,N) = 0.0
            CH_TILE(L,N) = 0.0
            Z0H_TILE(L,N) = 0.0
            Z0M_EFF_TILE(L,N) = 0.0
            RHOKPM(L,N) = 0.0
            RHOKPM_POT(L,N) = 0.0
            RADNET_TILE(L,N) = 0.0
            RHOKH_1(L,N) = 0.0
            FRACA(L,N) = 0.0
            RHOKM_1_TILE(L,N) = 0.0
            RESFT(L,N) = 0.0
            RESFS(L,N) = 0.0
          ENDIF
        ENDDO
      ENDDO

!  Set value of TINY_SNOW for > test(s)
!  Could use TINY (~10**-308) or EPSILON (~10**-16)
!  TINY turns out to be unsafe, but SQRT(TINY) (~10**-154) is ok
      TINY_SNOW = SQRT(TINY(SNOW_TILE(1,1)))
!!!      TINY_SNOW = EPSILON(SNOW_TILE(1,1))

!-----------------------------------------------------------------------
!!  1. Index array for sea-ice
!-----------------------------------------------------------------------

      NSICE = 0
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        IF ( ICE_FRACT(I,J) >  0.0 .AND. FLANDG(I,J) <  1.0 ) THEN
          NSICE = NSICE + 1
          SICE_INDEX(NSICE,1) = I
          SICE_INDEX(NSICE,2) = J
        ENDIF
       ENDDO
      ENDDO
!-----------------------------------------------------------------------
!!  1.1 Calculate height of lowest model level relative to sea.
!-----------------------------------------------------------------------

      DO J=1,ROWS
        DO I=1,ROW_LENGTH
          Z1_TQ_SEA(I,J)=Z1_TQ(I,J)
          IF(L_CTILE.AND.FLANDG(I,J) >  0.0.AND.FLANDG(I,J) <  1.0)     &
     &      Z1_TQ_SEA(I,J)=Z1_TQ(I,J)+Z_LAND(I,J)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!!  2.  Calculate QSAT values required later.
!-----------------------------------------------------------------------

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        RHOSTAR(I,J) = PSTAR(I,J) /                                     &
     &    ( R*(FLANDG(I,J)*TSTAR_LAND(I,J) +                            &
     &    (1.-FLANDG(I,J))*TSTAR_SSI(I,J)) )
!                        ... surface air density from ideal gas equation
       ENDDO
      ENDDO
! DEPENDS ON: qsat_mix
      CALL QSAT_mix(QS1,TL_1,PSTAR,ROW_LENGTH*ROWS,lq_mix_bl)
! DEPENDS ON: qsat_mix
      CALL QSAT_mix(QSTAR_SEA,TSTAR_SEA,PSTAR,ROW_LENGTH*ROWS           &
     & ,lq_mix_bl)
! DEPENDS ON: qsat_mix
      CALL QSAT_mix(QSTAR_ICE,TSTAR_SICE,PSTAR,ROW_LENGTH*ROWS          &
     & ,lq_mix_bl)
      DO L=1,LAND_PTS
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        PSTAR_LAND(L) = PSTAR(I,J)
      ENDDO
      DO N=1,NTILES
! DEPENDS ON: qsat_mix
        CALL QSAT_mix(QSTAR_TILE(1,N),TSTAR_TILE(1,N),                  &
     &            PSTAR_LAND,LAND_PTS,lq_mix_bl)
      ENDDO

!-----------------------------------------------------------------------
!!  3. Calculation of transfer coefficients and surface layer stability
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!!  3.1 Calculate neutral roughness lengths
!-----------------------------------------------------------------------

! Sea, sea-ice leads, sea-ice and marginal ice zone
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        Z0_MIZ(I,J) = Z0MIZ
        Z0_ICE(I,J) = Z0SICE
        RIB_SEA(I,J) = 0.
        RIB_ICE(I,J) = 0.
        DB_SEA(I,J) = 0.
        DB_ICE(I,J) = 0.
       ENDDO
      ENDDO
!
      IF (ISeaZ0T == SurfDivZ0T) THEN

!
!       Composite formulation for thermal roughness lengths,
!       incoporating the smooth aerodynamic limit for low
!       wind speeds and a value based on surface divergence
!       theory for higher wind speeds.
!
!       The friction velocity is diagnosed in the surface
!       transfer scheme, using z0m from the previous time-step.
!       z0[T,q] is also required but depends on u_* in this
!       scheme. For practical purposes, it is sufficient to
!       infer it from u_* determined from z0m, but in general
!       a given value of z0m corresponds to two values of
!       u_*, so we need to know whether we are on the low or
!       high wind speed side of the minimum value of z0m.
!       If we are on the high side, z0[T,q] will be inversely
!       proportional to z0m, but on the low side it may follow
!       this relationship, or be aerodynamically smooth. In
!       the smooth case we iterate u_* from the expression for
!       the roughness length and then take the maximum of the
!       smooth and high-wind expressions for z0[T,q]. An
!       iteration for the low-wind value of u_*, USTR_L is
!       carried out. This will converge to the correct limit only
!       on the low-wind side of the minimum, and the standard
!       criterion that the gradient of a fixed-point iteration
!       should be less than 1 in modulus gievs a more precise
!       boundary, TOL_USTR_L. For consistency with earlier versions
!       of the modset, hard-wired values are retained for the
!       operational value of Charnock's parameter. An additional
!       check is required, since z0m can be large at low or at
!       high wind-speeds. This is less precise and a fixed
!       value of 0.07 is used to test USTR_N, which was determined
!       by inspection of a graph of z0m against u_*: it is
!       unlikely that this will need to be changed unless
!       Charnock's constant is greatly altered.

!	YMA: compute TOL_USTR_L with the new Charnock expression
          TOL_USTR_L = 0.75*(1.54E-6*G/(2.0*CK4Z0ML))**0.33333
          TOL_USTR_N = 0.07
!
        DO j=1, rows
          DO I=1, row_length
!           We need to infer u_* from Z0M.
            IF (VSHR_SSI(I,j)  > 0.0) THEN

!	if (i==6 .and. j==1) then
!	  print *, 'dhb1: VSHR_SSI(I,j)= ', VSHR_SSI(I,j)  
!	  print *, 'dhb1: RT4CK,UMD4CK= ',RT4CK,UMD4CK
!	endif

!	YMA: set new Charnock expression to be used in new z0h form
         XCK4Z0M=RT4CK*(VSHR_SSI(I,J)-UMD4CK)
         EXPOSI=EXP(XCK4Z0M)
         EXNEGA=1./EXPOSI   !EXP(-XCK4Z0M)

!        if (i==6 .and. j==1) then
!          print *, 'dhb2: XCK4Z0M= ', XCK4Z0M
!        endif

         CK4Z0M=((EXPOSI-EXNEGA)/(EXPOSI+EXNEGA))*0.5+0.5

!        if (i==6 .and. j==1) then
!          print *, 'dhb3: XCK4Z0M= ', CK4Z0M
!          print *, 'dhb3: CK4Z0MH,CK4Z0ML= ', CK4Z0MH,CK4Z0ML
!        endif

         CK4Z0M=(CK4Z0MH-CK4Z0ML)*CK4Z0M + CK4Z0ML

!        if (i==6 .and. j==1) then
!          print *, 'dhb4: CK4Z0M= ', CK4Z0M
!        endif

!        if (i==6 .and. j==1) then
!          print *, 'dhb5: VKMAN,VSHR_SSI= ', VKMAN,VSHR_SSI(I,j)
!          print *, 'dhb5: Z1_UV,Z0MSEA  = ', Z1_UV(I,j),Z0MSEA(I,j) 
!        endif

!             Compute u_* using neutral stratification.
!             stratification.
              USTR_N = VKMAN * VSHR_SSI(I,j) /                        &
     &          LOG(Z1_UV(I,j) / Z0MSEA(I,j) )
!             Compute u_* using low wind approximation.
              USTR_L = 1.54E-06 /  Z0MSEA(I,j) - 1.0E-05
!             Since Z0M could be large for low and high u_*, we use
!             USTR_N as an additional check on the regime.

!dhb: filtering out too small USTR_L, esp assuring NO neg value ---------------------
              !if (USTR_L < 1.e-6) then
              if (USTR_L < 1.e-8) then
                write(*,'(a,2i5,e18.8)') 'dhb01: Warning--I,j, USTR_L= ', I,j, USTR_L
                !USTR_L = max(USTR_L, 1.e-6)
                USTR_L = max(USTR_L, 1.e-8)
              endif
!------------------------------------------------------------------------------------

!        if (i==6 .and. j==1) then
!          print *, 'dhb6: USTR_N, USTR_L= ', USTR_N, USTR_L
!          print *, 'dhb6: TOL_USTR_N, TOL_USTR_L= ', TOL_USTR_N, TOL_USTR_L
!        endif

              IF ( (USTR_N < TOL_USTR_N) .AND.                        &
     &             (USTR_L < TOL_USTR_L) ) THEN
!               Iterate u_* for low winds.
                DO JITS=1, 5

!        if (i==6 .and. j==1) then
!          print *, 'dhb7xxx: Itererate No JITS= ',JITS
!          print *, 'dhb7xxx: Z0MSEA,CK4Z0M= ',Z0MSEA(I,j),CK4Z0M
!          print *, 'dhb7xxx: G= ',G
!          print *, 'dhb7xxx: BF: USTR_L= ',USTR_L
!        endif

!		YMA: the new Charnock parameter is used
                  USTR_L=1.54E-06/(Z0MSEA(I,j)-(CK4Z0M/G)*USTR_L**2)  &
     &              -1.0E-05

!dhb: filtering out too small USTR_L, esp assuring NO neg value ---------------------
              !if (USTR_L < 1.e-6) then
              if (USTR_L < 1.e-8) then
                write(*,'(a,3i5,e18.8)') 'dhb02: I,j,JITS,USTR_L= ', I,j,JITS,USTR_L
                !USTR_L = max(USTR_L, 1.e-6)
                USTR_L = max(USTR_L, 1.e-8)
                exit
              endif
!------------------------------------------------------------------------------------

!        if (i==6 .and. j==1) then
!          print *, 'dhb7xxx: AF: USTR_L= ',USTR_L
!        endif
                ENDDO
!               Take the maximum of the smooth and high-wind values.
!               A lower limit is imposed on the friction velocity to
!               allow for the exceptional case of very low winds: the
!               value of 10^-5 is the same as the limit for the momentum
!               roughness length.
              ELSE
!               Take the high-wind value, but limit it to the molecular
!               mean free path (we should not hit this limit
!               in practice).

!        if (i==6 .and. j==1) then
!          print *, 'dhb8: Z0MSEA,CK4Z0M,G== ',Z0MSEA(I,j),CK4Z0M,G
!        endif

!		YMA: the new Charnock parameter is used
                USTR_L=sqrt(Z0MSEA(I,j)/(CK4Z0M/G))

!        if (i==6 .and. j==1) then
!          print *, 'dhb9: USTR_L= ',USTR_L
!        endif

              ENDIF

!	YMA: set new z0h from
!        if (i==6 .and. j==1) then
!          print *, 'dhb10: Z0MSEA,USTR_L= ',Z0MSEA(I,j),USTR_L
!        endif
        RT4Z0H = Z0MSEA(I,j)*USTR_L/1.40E-5
!        if (i==6 .and. j==1) then
!          print *, 'dhb11XXX: RT4Z0H= ',RT4Z0H
!        endif
        ALGRRDIF4Z0H  = 1.5 * (ALOG10(RT4Z0H) - 0.25 )
!        if (i==6 .and. j==1) then
!          print *, 'dhb12XXX: ALGRRDIF4Z0H= ',ALGRRDIF4Z0H
!        endif
        EXPOSI=EXP(ALGRRDIF4Z0H)
        EXNEGA=1./EXPOSI !EXP(-ALGRRDIF4Z0H)
        SH4Z0H = (EXPOSI-EXNEGA)/(EXPOSI+EXNEGA)
!        if (i==6 .and. j==1) then
!          print *, 'dhb13XXX: SH4Z0H= ',SH4Z0H
!        endif
        ALGZ0T4Z0H = (-3.8 - (-4.8))*(-SH4Z0H*0.5+0.5)+(-4.8)
!        if (i==6 .and. j==1) then
!          print *, 'dhb14XXX: RT4Z0H,ALGZ0T4Z0H= ',RT4Z0H,ALGZ0T4Z0H
!        endif
        ALGZ0T4Z0H = ALGZ0T4Z0H + ((ALOG10(RT4Z0H) - (-1.))/(3.-(-1.)))*(-0.20)
!        if (i==6 .and. j==1) then
!          print *, 'dhb15XXX: ALGZ0T4Z0H= ',ALGZ0T4Z0H
!        endif

        Z0H_SEA(I,J) = MIN(1.1E-4, 10.**ALGZ0T4Z0H)

!        if (i==6 .and. j==1) then
!          print *, 'dhb16XXX: Z0H_SEA= ',Z0H_SEA(I,J)
!        endif

!	if((i==6).and.(j==1)) then
!	print*,"WRONGYMA:i,j,Z0H_SEA(I,J)=", i,j,Z0H_SEA(I,J)
!	print*,"Z0MSEA(I,j), USTR_L=", Z0MSEA(I,j), USTR_L
!	endif

            ENDIF
          ENDDO
        ENDDO
!
      ELSE IF (ISeaZ0T == Fixed_Z0T) THEN
!
!       Use a fixed thermal roughness length.
        Z0H_SEA(1:row_length,1:rows) = Z0HSEA
!
      ENDIF
!

      IF ( L_spec_z0 ) THEN
         DO J = 1, rows
         DO I = 1, row_length
           IF (Z0H_SCM(I,J)  >   0.0) THEN

             ! Set Z0H from SCM namelist
             ! (if specified) for sea points
             Z0H_SEA(I,J) = Z0H_SCM(I,J)

           ENDIF  ! Z0H_SCM
         ENDDO  ! I
         ENDDO  ! J
      ENDIF

! Land tiles
! Z0_TILE contains the appropriate value for land-ice points, but has to
! be modified for snow-cover on non-land-ice points
      DO N=1,NTILES
!CDIR NODEP
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
!!! Maybe this ought to be > TINY_SNOW too, but no problems seen
          IF ( SNOW_TILE(L,N) >  0. .AND. SMVCST(L) /= 0. ) THEN
            Z0 = Z0_TILE(L,N) - 4.0E-4*SNOW_TILE(L,N)
            ZETA1 = MIN( 5.0E-4 , Z0_TILE(L,N)  )
            Z0M_TILE(L,N) = MAX( ZETA1 , Z0 )
          ELSE
            Z0M_TILE(L,N) = Z0_TILE(L,N)
          ENDIF
          Z0H_TILE(L,N) = Z0H_Z0M(N)*Z0M_TILE(L,N)
          RIB_TILE(L,N) = 0.
          DB_TILE(L,N) = 0.
        ENDDO
      ENDDO

      DO N=1,NTILES
! DEPENDS ON: sf_orog
        CALL SF_OROG (                                                  &
     &   ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS(N),                          &
     &   LAND_INDEX,TILE_INDEX(1,N),                                    &
     &   FORMDRAG,FD_stab_dep,OROG_DRAG_PARAM,LTIMER,                   &
     &   HO2R2_OROG,RIB_TILE(1,N),SIL_OROG_LAND,Z0M_TILE(1,N),Z1_UV,    &
     &   WIND_PROFILE_FACTOR(1,N),Z0M_EFF_TILE(1,N)                     &
     &   )
      ENDDO

! EAK
      conv_rain = 0. ! convective rain and snow calculated after 
                     ! call to ni_bl_ctl
      conv_snow = 0.

!      print *,'sf_exch bef cable_Exum l_cable',l_cable
!EAK  SW_DOWN calculations for checking only !!!
      SW_DOWN = 0.
      SW_DOWN = (surf_down_sw(:,:,1)+surf_down_sw(:,:,2) + &
         surf_down_sw(:,:,3)+surf_down_sw(:,:,4))*cos_zenith_angle(:,:)

      if( l_cable ) then
! DEPENDS ON: cable_explicit_driver.o
         CALL cable_explicit_driver(                                        &
                   row_length, rows, land_pts, ntiles,    &
                   npft, sm_levels, timestep, lat, long, land_index,       &
                   tile_frac, tile_pts, tile_index,&
                   bexp, hcon, satcon, sathh,       &
                   smvcst, smvcwt, smvccl, albsoil, snow_tile, snow_rho1l,           &
                   snage_tile, isnow_flg3l, snow_rho3l, snow_cond, snow_depth3l,     &
                   snow_tmp3l, snow_mass3l, sw_down, lw_down, cos_zenith_angle,      &
                   surf_down_sw, ls_rain, ls_snow, tl_1, qw_1, vshr_land, pstar,     &
                   z1_tq, z1_uv, rho_water, L_tile_pts, canopy, Fland,          &
! rml 2/7/13 pass 3d co2 through to cable if required
                   CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,         &
                   sthu_tile, smcl_tile, sthf_tile, sthu, tsoil_tile,                &
                   canht_ft, lai_ft, sin_theta_latitude, dzsoil, &
                   LAND_MASK, FTL_TILE_CAB, FTL_CAB,    &
                   FTL_TILE, FQW_TILE, LE_TILE_CAB, LE_CAB, TSTAR_TILE,              &
                   TSTAR_TILE_CAB, TSTAR_CAB, U_S, U_S_STD_TILE, U_S_CAB, CH_CAB,    &
                   CD_CAB, CD_TILE, CH_TILE, RADNET_TILE, FRACA, rESFS, RESFT, &
                   Z0H_TILE, Z0M_TILE, RECIP_L_MO_TILE, EPOT_TILE,  &
                   ! Lestevens Sept 2012 - CasaCNP variables
                   CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                   &
                   NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                         &
     ! Lestevens 23apr13
                   NPP_FT_ACC,RESP_W_FT_ACC,                                      &
                   ! end step of experiment, this step, step width, processor num
                   endstep, timestep_number, mype, Lai_Ma &
                 & )
      endif

! EAK
      IF( l_cable ) THEN
        CD_STD = CD_TILE
        V_S_TILE = U_S_STD_TILE
        V_S_STD = U_S_STD_TILE
      ENDIF ! l_cable

!-----------------------------------------------------------------------
! Calculate RESFT with neutral CH and EPDT=0 for use in calculation
! of Richardson number. RESFT=1 for snow and land-ice.
!-----------------------------------------------------------------------

      DO N=1,NTILES
!CDIR NODEP
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          ZETAM = LOG ( (Z1_UV(I,J) + Z0M_TILE(L,N))/max(1.e-7,Z0M_TILE(L,N)) )
          ZETAH = LOG ( (Z1_TQ(I,J) + Z0M_TILE(L,N))/max(1.e-7,Z0H_TILE(L,N)) )
          CHN(L) = (VKMAN/ZETAH)*(VKMAN/ZETAM)*WIND_PROFILE_FACTOR(L,N)
          DQ(L) = QW_1(I,J) - QSTAR_TILE(L,N)
          EPDT(L) = 0.0
!EAK
          if (l_cable) RESFT(L,N) =  min(1.,FLAKE(L,N) + (1. - FLAKE(L,N)) *        &
     &                ( FRACA(L,N) + (1. - FRACA(L,N))*RESFS(L,N) )) 
        ENDDO
! DEPENDS ON: sf_resist
      IF( .NOT.  l_cable ) THEN

        CALL SF_RESIST (                                                &
     &   ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS(N),                          &
     &   LAND_INDEX,TILE_INDEX(1,N),                                    &
     &   CANOPY(1,N),CATCH(1,N),CHN,DQ,EPDT,FLAKE(1,N),GC(1,N),         &
     &   SNOW_TILE(1,N),VSHR_LAND,FRACA(1,N),RESFS(1,N),RESFT(1,N),     &
     &   LTIMER                                                         &
     &   )
      ENDIF ! .NOT.  l_cable 
      ENDDO

! RESFT < 1 for snow on needleleaf tile if canopy snow model used
      IF( .NOT.  l_cable ) THEN
        IF (NTILES >  1 .AND. CAN_MODEL == 4) THEN
        N = 2   ! NLT tile
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          IF (SNOW_TILE(L,N)  > TINY_SNOW) THEN
            CE_D_SH(L) = 0.02*(CATCH_SNOW(L)/SNOW_TILE(L,N))**0.4       &
     &                   * 2.06E-5*(TM/TSTAR_TILE(L,N))**1.75           &
     &                   * (1.79 + 3*SQRT(VSHR_LAND(I,J)))
            GC(L,N) = 3*SNOW_TILE(L,N)*CE_D_SH(L)/(2*RHO_ICE*5E-4**2)
            FRACA(L,N) = 0.
            RESFS(L,N) = GC(L,N)/(GC(L,N) + CHN(L)*VSHR_LAND(I,J))
            RESFT(L,N) = RESFS(L,N)
          ENDIF
        ENDDO
       ENDIF
      ENDIF ! .NOT.  l_cable 


!-----------------------------------------------------------------------
!!  3.2 Calculate bulk Richardson number for the lowest model level.
!-----------------------------------------------------------------------

! Sea, sea-ice and sea-ice leads
! DEPENDS ON: sf_rib_sea
      CALL SF_RIB_SEA (                                                 &
     & ROW_LENGTH,ROWS,FLANDG,NSICE,SICE_INDEX,                         &
     & BQ_1,BT_1,ICE_FRACT,QSTAR_ICE,QSTAR_SEA,QW_1,TL_1,TSTAR_SICE,    &
     & TSTAR_SEA,VSHR_SSI,Z0_ICE,Z0H_SEA,Z0_ICE,Z0MSEA,Z1_TQ_SEA,Z1_UV, &
     & RIB_SEA,RIB_ICE,DB_SEA,DB_ICE,LTIMER                             &
     & )

! Land tiles
      DO N=1,NTILES
! DEPENDS ON: sf_rib_land
        CALL SF_RIB_LAND (                                              &
     &   ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS(N),                          &
     &   LAND_INDEX,TILE_INDEX(1,N),                                    &
     &   BQ_1,BT_1,QSTAR_TILE(1,N),QW_1,RESFT(1,N),TL_1,                &
     &   TSTAR_TILE(1,N),VSHR_LAND,Z0H_TILE(1,N),Z0M_TILE(1,N),         &
     &   Z1_TQ,Z1_UV,                                                   &
     &   RIB_TILE(1,N),DB_TILE(1,N),LTIMER                              &
     &   )
      ENDDO

!-----------------------------------------------------------------------
!!  3.3 Calculate stability corrected effective roughness length.
!!  Stability correction only applies to land points.
!-----------------------------------------------------------------------

      DO N=1,NTILES
! DEPENDS ON: sf_orog
        CALL SF_OROG (                                                  &
     &   ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS(N),                          &
     &   LAND_INDEX,TILE_INDEX(1,N),                                    &
     &   FORMDRAG,FD_stab_dep,OROG_DRAG_PARAM,LTIMER,                   &
     &   HO2R2_OROG,RIB_TILE(1,N),SIL_OROG_LAND,Z0M_TILE(1,N),Z1_UV,    &
     &   WIND_PROFILE_FACTOR(1,N),Z0M_EFF_TILE(1,N)                     &
     &   )
      ENDDO

!-----------------------------------------------------------------------
!!  3.4 Calculate CD, CH via routine FCDCH.
!!      Note that these are returned as the dimensionless surface
!!      exchange coefficients.
!-----------------------------------------------------------------------

! Sea-ice
! DEPENDS ON: fcdch_sea
      CALL FCDCH_SEA(ROW_LENGTH,ROWS,COR_MO_ITER,FLANDG,                &
     &               DB_ICE,VSHR_SSI,Z0_ICE,Z0_ICE,                     &
     &               ZH,Z1_UV,Z1_TQ_SEA,                                &
     &               CD_ICE,CH_ICE,V_S_ICE,                             &
     &               RECIP_L_MO_ICE,LTIMER)

! Marginal Ice Zone
! DEPENDS ON: fcdch_sea
      CALL FCDCH_SEA(ROW_LENGTH,ROWS,COR_MO_ITER,FLANDG,                &
     &               DB_ICE,VSHR_SSI,Z0_MIZ,Z0_MIZ,                     &
     &               ZH,Z1_UV,Z1_TQ_SEA,                                &
     &               CD_MIZ,CH_MIZ,V_S_MIZ,                             &
     &               RECIP_L_MO_MIZ,LTIMER)

! Sea and sea-ice leads
! DEPENDS ON: fcdch_sea
      CALL FCDCH_SEA(ROW_LENGTH,ROWS,COR_MO_ITER,FLANDG,                &
     &               DB_SEA,VSHR_SSI,Z0MSEA,Z0H_SEA,                    &
     &               ZH,Z1_UV,Z1_TQ_SEA,                                &
     &               CD_SEA,CH_SEA,V_S_SEA,                             &
     &               RECIP_L_MO_SEA,LTIMER)

! Land tiles
      IF( .NOT.  l_cable ) THEN
       DO N=1,NTILES
! DEPENDS ON: fcdch_land
        CALL FCDCH_LAND (                                               &
     &   ROW_LENGTH,ROWS,COR_MO_ITER,LAND_PTS,TILE_PTS(N),              &
     &   TILE_INDEX(1,N),LAND_INDEX,                                    &
     &   DB_TILE(1,N),VSHR_LAND,                                        &
     &   Z0M_EFF_TILE(1,N),Z0H_TILE(1,N),ZH,                            &
     &   Z1_UV,Z1_TQ,WIND_PROFILE_FACTOR(1,N),                          &
     &   CD_TILE(1,N),CH_TILE(1,N),CD_STD(1,N),                         &
     &   V_S_TILE(1,N),V_S_STD(1,N),RECIP_L_MO_TILE(1,N),               &
     &   U_S_ITER_TILE(1,N),L_DUST,                                     &
     &   LTIMER                                                         &
     &   )
      ENDDO

      IF ( COR_MO_ITER == Use_Correct_Ustar ) THEN
!       Use correct "standard" ustar
        U_S_STD_TILE(:,:) = V_S_STD(:,:)
      ELSE
!       Use ustar from mid-iteration
        U_S_STD_TILE(:,:) = U_S_ITER_TILE(:,:)
      END IF

! MRD - is the end if in the correct place or should it be after the enddo above?
      ENDIF ! .NOT.  l_cable

!-----------------------------------------------------------------------
!!  4.1 Recalculate RESFT using "true" CH and EPDT for land tiles
!-----------------------------------------------------------------------

      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          DQ(L) = QW_1(I,J) - QSTAR_TILE(L,N)
          EPDT(L) = - RHOSTAR(I,J)*CH_TILE(L,N)*VSHR_LAND(I,J)          &
     &      *DQ(L)*TIMESTEP
        ENDDO
! DEPENDS ON: sf_resist
      IF( .NOT.  l_cable ) THEN
        CALL SF_RESIST (                                                &
     &   ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS(N),                          &
     &   LAND_INDEX,TILE_INDEX(1,N),                                    &
     &   CANOPY(1,N),CATCH(1,N),CH_TILE(1,N),DQ,EPDT,FLAKE(1,N),        &
     &   GC(1,N),SNOW_TILE(1,N),VSHR_LAND,FRACA(1,N),                   &
     &   RESFS(1,N),RESFT(1,N),                                         &
     &   LTIMER)
      ENDIF !  IF( .NOT.  l_cable 
      ENDDO

      IF( .NOT.  l_cable ) THEN
       IF (NTILES >  1 .AND. CAN_MODEL == 4) THEN
        N = 2
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
!!! Maybe this ought to be > TINY_SNOW too, but no problems seen
          IF (SNOW_TILE(L,N)  >   0.) THEN
            FRACA(L,N) = 0.
            RESFS(L,N) = GC(L,N) /                                      &
     &                  (GC(L,N) + CH_TILE(L,N)*VSHR_LAND(I,J))
            RESFT(L,N) = RESFS(L,N)
          ENDIF
        ENDDO
       ENDIF
      ENDIF !  IF( .NOT.  l_cable 

!-----------------------------------------------------------------------
! Calculate gridbox-means of transfer coefficients.
!-----------------------------------------------------------------------

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        CD_LAND(I,J) = 0.
        CD_SSI(I,J) = 0.
        CH_SSI(I,J) = 0.
        CD(I,J) = 0.
        CH(I,J) = 0.
        IF (L_DUST) CD_STD_DUST(I,J) = 0.

! Sea and sea-ice
        IF (FLANDG(I,J) <  1.0 ) THEN
          IF ( ICE_FRACT(I,J)  <   0.7 ) THEN
            CD_SSI(I,J) = ( ICE_FRACT(I,J)*CD_MIZ(I,J) +                &
     &              (0.7-ICE_FRACT(I,J))*CD_SEA(I,J) ) / 0.7  ! P2430.5
            CH_SSI(I,J) = ( ICE_FRACT(I,J)*CH_MIZ(I,J) +                &
     &              (0.7-ICE_FRACT(I,J))*CH_SEA(I,J) ) / 0.7  ! P2430.4
          ELSE
            CD_SSI(I,J) = ( (1.0-ICE_FRACT(I,J))*CD_MIZ(I,J) +          &
     &              (ICE_FRACT(I,J)-0.7)*CD_ICE(I,J) ) / 0.3  ! P2430.7
            CH_SSI(I,J) = ( (1.0-ICE_FRACT(I,J))*CH_MIZ(I,J) +          &
     &              (ICE_FRACT(I,J)-0.7)*CH_ICE(I,J) ) / 0.3  ! P2430.7
          ENDIF
          CD(I,J)=(1.-FLANDG(I,J))*CD_SSI(I,J)
          CH(I,J)=(1.-FLANDG(I,J))*CH_SSI(I,J)
          IF (L_DUST) CD_STD_DUST(I,J)=(1.-FLANDG(I,J))*CD_SSI(I,J)
        ENDIF

       ENDDO
      ENDDO

! Land tiles
      DO N=1,NTILES
!CDIR NODEP
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          CD_LAND(I,J) = CD_LAND(I,J) + TILE_FRAC(L,N)*CD_TILE(L,N)
          CD(I,J) = CD(I,J) + FLANDG(I,J)*TILE_FRAC(L,N)*CD_TILE(L,N)
          CH(I,J) = CH(I,J) + FLANDG(I,J)*TILE_FRAC(L,N)*CH_TILE(L,N)
          IF (L_DUST) CD_STD_DUST(I,J) = CD_STD_DUST(I,J) +             &
     &                  FLANDG(I,J)*TILE_FRAC(L,N)*CD_STD(L,N)

        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!!  4.3 Calculate the surface exchange coefficients RHOK(*) and
!       resistances for use in Sulphur Cycle
!       (Note that CD_STD, CH and VSHR should never = 0)
!     RHOSTAR * CD * VSHR stored for diagnostic output before
!     horizontal interpolation.
!-----------------------------------------------------------------------

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        RHO_ARESIST(I,J) = 0.
         RHO_ARESIST_LAND(I,J)=0.0
        ARESIST(I,J) = 0.
        RESIST_B(I,J) = 0.
        RHOKM_LAND(I,J) = 0.
        RHOKM_SSI(I,J) = 0.

! Sea and sea-ice
        IF ( FLANDG(I,J) <  1.0 ) THEN
          RHOKM_SSI(I,J) = RHOSTAR(I,J)*CD_SSI(I,J)*VSHR_SSI(I,J)
!                                                          ! P243.124
          RHOKH_1_SICE(I,J) = RHOSTAR(I,J)*CH_SSI(I,J)*VSHR_SSI(I,J)
!                                                           ! P243.125
          RHO_ARESIST(I,J) = RHOSTAR(I,J)*CD_SSI(I,J)*VSHR_SSI(I,J)
          ARESIST(I,J) =  1. / (CD_SSI(I,J) * VSHR_SSI(I,J))
          RESIST_B(I,J)= (CD_SSI(I,J)/CH_SSI(I,J) - 1.0) * ARESIST(I,J)
        ENDIF

       ENDDO
      ENDDO

! Land tiles
      DO N=1,NTILES

        DO L=1,LAND_PTS
          RHO_ARESIST_TILE(L,N) = 0.
          ARESIST_TILE(L,N) = 0.
          RESIST_B_TILE(L,N) = 0.
        ENDDO
!CDIR NODEP
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          RHOKM_1_TILE(L,N) = RHOSTAR(I,J)*CD_TILE(L,N)*VSHR_LAND(I,J)
!                                                         ! P243.124
          RHOKM_LAND(I,J) = RHOKM_LAND(I,J) +                           &
     &           TILE_FRAC(L,N)*RHOKM_1_TILE(L,N)
          RHOKH_1(L,N) = RHOSTAR(I,J)*CH_TILE(L,N)*VSHR_LAND(I,J)
!                                                         ! P243.125
          RHO_ARESIST_TILE(L,N) = RHOSTAR(I,J) * CD_STD(L,N)            &
     &                * VSHR_LAND(I,J)
          ARESIST_TILE(L,N) = 1. / ( CD_STD(L,N) * VSHR_LAND(I,J) )
          RESIST_B_TILE(L,N) = ( CD_STD(L,N)/CH_TILE(L,N) - 1.0 ) *     &
     &                                                 ARESIST_TILE(L,N)
          IF (RESIST_B_TILE(L,N)  <   0.) RESIST_B_TILE(L,N) = 0.
          RHO_ARESIST_LAND(I,J) = RHO_ARESIST_LAND(I,J) +               &
     &                     TILE_FRAC(L,N)*RHO_ARESIST_TILE(L,N)
        ENDDO
      ENDDO
!CDIR NODEP
      DO L=1,LAND_PTS
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        RHO_ARESIST(I,J) = FLANDG(I,J)*RHO_ARESIST_LAND(I,J) +          &
     &                      (1.0-FLANDG(I,J))*RHO_ARESIST(I,J)
        ARESIST(I,J) = RHOSTAR(I,J) / RHO_ARESIST(I,J)
      ENDDO

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        RHOKM_1(I,J)= FLANDG(I,J) * RHOKM_LAND(I,J) +                   &
     &                     (1.0-FLANDG(I,J)) * RHOKM_SSI(I,J)
        RHO_CD_MODV1(I,J) = RHOKM_1(I,J)  ! diagnostic required for VAR
       ENDDO
      ENDDO


!-----------------------------------------------------------------------
!!  Calculate surface layer resistance for mineral dust
!-----------------------------------------------------------------------
       IF (L_DUST) THEN

         DO J = 1,ROWS
           DO I = 1,ROW_LENGTH
             DO IDIV=1,NDIV
               R_B_DUST(I,J,IDIV)=0.
             ENDDO !IDIV
             VSHR(I,J)=(1.0 - FLANDG(I,J)) * VSHR_SSI(I,J) +            &
     &              FLANDG(I,J) * VSHR_LAND(I,J)
           ENDDO !I
         ENDDO !J

! DEPENDS ON: dustresb
         CALL DUSTRESB (                                                &
     &  ROW_LENGTH,ROWS,                                                &
     &  PSTAR,TSTAR,RHOSTAR,ARESIST,VSHR,CD_STD_DUST,                   &
     &  R_B_DUST                                                        &
     &  )

       ENDIF !(L_DUST)
!-----------------------------------------------------------------------
!!  Calculate local and gridbox-average surface fluxes of heat and
!!  moisture.
!-----------------------------------------------------------------------

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        FTL_1(I,J) = 0.
        FQW_1(I,J) = 0.
       ENDDO
      ENDDO

! Sea and sea-ice
! DEPENDS ON: sf_flux_sea
      CALL SF_FLUX_SEA (                                                &
     & ROW_LENGTH,ROWS,NSICE,SICE_INDEX,FLANDG,                         &
     & ICE_FRACT,QS1,QSTAR_ICE,QSTAR_SEA,QW_1,RADNET_SICE,RHOKH_1_SICE, &
     & TI,TL_1,TSTAR_SICE,TSTAR_SEA,Z0_ICE,Z0_ICE,Z0H_SEA,Z0MSEA,       &
     & Z1_TQ_SEA,SeaSalinityFactor,                                     &
     & ALPHA1_SICE,ASHTF,E_SEA,FQW_ICE,FQW_1,FTL_ICE,FTL_1,H_SEA,       &
     & RHOKPM_SICE,LTIMER                                               &
     & )

! Land tiles

! EAK
      IF(  l_cable ) THEN
       DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          D_T = TSTAR_TILE(L,N) - TL_1(I,J)
          IF (D_T  >   0.05 .OR. D_T  <   -0.05) THEN
            ALPHA1(L,N) = (QSTAR_TILE(L,N) - QS1(I,J)) / D_T
          ELSEIF (TL_1(I,J)  >   TM) THEN
            ALPHA1(L,N) = EPSILON*LC*QS1(I,J)*                              &
     &        (1. + C_VIRTUAL*QS1(I,J)) / ( R*TL_1(I,J)*TL_1(I,J))
          ELSE
            ALPHA1(L,N) = EPSILON*LS*QS1(I,J)*                              &
     &        (1. + C_VIRTUAL*QS1(I,J)) / ( R*TL_1(I,J)*TL_1(I,J))
          ENDIF
! This needs to be look at:
          ASHTF_TILE(L,N) = 2.0 * HCONS(L) / DZSOIL
          IF (SNOW_TILE(L,N) >  0.0 .AND. SMVCST(L) /= 0.) THEN
            DS_RATIO = 2.0 * SNOW_TILE(L,N) / (RHO_SNOW * DZSOIL)
            IF (DS_RATIO <= 1.0) THEN
              ASHTF_TILE(L,N) =  ASHTF_TILE(L,N) /                         &
     &                    (1. + DS_RATIO*(HCONS(L)/SNOW_HCON - 1.))
            ELSE
              ASHTF_TILE(L,N) =  ASHTF_TILE(L,N)*SNOW_HCON / HCONS(L)
            ENDIF
          ENDIF
!          ASHTF_TILE(L,N) = (1.-VFRAC_TILE(L,N))*ASHTF_TILE(L,N)+          &
!     &           CANHC_TILE(L,N)      &
!     &          /TIMESTEP + 4*(1. + VFRAC_TILE(L,N))*SBCON*TS1(L)**3
          LH = LC
          IF (SNOW_TILE(L,N)  >   0.) LH = LS
          RHOKPM(L,N) = RHOKH_1(L,N) / ( ASHTF_TILE(L,N)  +                &
     &                  RHOKH_1(L,N)*(LH*ALPHA1(L,N)*RESFT(L,N) + CP))
          RHOKPM_POT(L,N)=RHOKH_1(L,N) / ( ASHTF_TILE(L,N)  +              &
     &                     RHOKH_1(L,N)*(LH*ALPHA1(L,N) + CP) )
!          RAD_REDUC = RADNET_TILE(L,N)-ASHTF_TILE(L,N)*(TL_1(I,J)-TS1(L)   &
!     &       + GRCP*(Z1_TQ(I,J) + Z0M_EFF_TILE(L,N) - Z0H_TILE(L,N)) )    &
!     &       + CANHC_TILE(L,N)*(TSTAR_TILE(L,N) - TS1(L)) / TIMESTEP
!          DQ1 = QS1(I,J) - QW_1(I,J) + GRCP*ALPHA1(L,N)*                  &
!     &          (Z1_TQ(I,J) + Z0M_EFF_TILE(L,N) - Z0H_TILE(L,N))
!          FQW_TILE(L,N) = RESFT(L,N)*RHOKPM(L,N)*(ALPHA1(L,N)*RAD_REDUC   &
!     &               + (CP*RHOKH_1(L,N) + ASHTF_TILE(L,N))*DQ1 )
!          EPOT_TILE(L,N) = RHOKPM_POT(L,N)*( ALPHA1(L,N)*RAD_REDUC          &
!     &               + (CP*RHOKH_1(L,N) + ASHTF_TILE(L,N))*DQ1 )
!          FTL_TILE(L,N) = RHOKPM(L,N)*(RAD_REDUC - LH*RESFT(L,N)*           &
!     &                  RHOKH_1(L,N)*DQ1)
          FTL_1(I,J) = FTL_1(I,J)+FLAND(L)*TILE_FRAC(L,N)*FTL_TILE(L,N)
          FQW_1(I,J) = FQW_1(I,J)+FLAND(L)*TILE_FRAC(L,N)*FQW_TILE(L,N)

        ENDDO
       ENDDO
!      print *,'in SFEXCH USECABLE',FQW_TILE,FTL_TILE
!      print *,'in SFEXCH USECABLE 1',FQW_1,FTL_1

      ELSE


      DO N=1,NTILES
        DO L = 1,LAND_PTS
          FTL_TILE(L,N) = 0.
          FQW_TILE(L,N) = 0.
          EPOT_TILE(L,N) = 0.
        ENDDO
      ENDDO

      DO N=1,NTILES
! DEPENDS ON: sf_flux_land
        CALL SF_FLUX_LAND (                                             &
     &   ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS(N),FLAND,                    &
     &   LAND_INDEX,TILE_INDEX(1,N),                                    &
     &   CANHC_TILE(1,N),DZSOIL,HCONS,QS1,QSTAR_TILE(1,N),QW_1,         &
     &   RADNET_TILE(1,N),RESFT(1,N),RHOKH_1(1,N),SMVCST,SNOW_TILE(1,N),&
     &   TILE_FRAC(1,N),TIMESTEP,TL_1,TS1,TSTAR_TILE(1,N),              &
     &   VFRAC_TILE(1,N),Z0H_TILE(1,N),Z0M_EFF_TILE(1,N),Z1_TQ,         &
     &   FQW_1,FTL_1,                                                   &
     &   ALPHA1(1,N),ASHTF_TILE(1,N),FQW_TILE(1,N),EPOT_TILE(1,N),      &
     &   FTL_TILE(1,N),                                                 &
     &   RHOKPM(1,N),RHOKPM_POT(1,N),LTIMER,ANTHROP_HEAT(N)             &
     & )
      ENDDO

       
      ENDIF    ! (  l_cable ) THEN

      IF(ntest>0) THEN

      print 184,ftl_1(4,2),fqw_1(4,2), &
                  ftl_tile(100,9),fqw_tile(100,9),tstar_tile(100,9)
      print 185,ftl_1(28,57),fqw_1(28,57), &
      ftl_tile(2132,1:5),ftl_tile(2132,8),fqw_tile(2132,1:5), &
      fqw_tile(2132,8),tstar_tile(2132,1:5),tstar_tile(2132,8)
      print 186,ftl_1(35,56),fqw_1(35,56), &
      ftl_tile(2080,1:5),ftl_tile(2080,8),fqw_tile(2080,1:5), &
      fqw_tile(2080,8),tstar_tile(2080,1:5),tstar_tile(2080,8)
184   format('sfexaflan',20f6.1)
185   format('sfexaflan',20f6.1)
186   format('sfexaflan',20f6.1)


!     print all out variables from sf_flux_land for 3 selected points
!      print 77,RHOKPM(100,:),RADNET_TILE(100,:),FQW_1(I,J),  &
!     &FQW_TILE(100,:),FTL_1(I,J),FTL_TILE(100,:),FLAND(L),TILE_FRAC(L,:)
!77    format(1x,'sfexch100',9f6.1,1x,9f5.0,1x,f6.1,1x,9f6.1,1x,f6.1,1x,9f6.1, &
!     & 1x,f5.2,9f5.2)
!     l=100 has i=4,j=2,N=9
!     l=1000 has i=13,j=29
!     l=1500 has i=94,j=44
!     l=100 has i=4,j=2
      print *,'l=100,i=4,j=2'
      l=100
      i=4
      j=2
      N=9
      print 21,l,RHOKPM(l,N),RHOKPM_POT(l,N),RADNET_TILE(l,N),   &
                 SW_DOWN(i,j),SW_TILE(l,N),LW_DOWN(i,j)
      print 22,l,RHOKH_1(l,N),SNOW_TILE(l,N),TL_1(i,j),TS1(l),QS1(i,j)
      print 23,l,FQW_1(i,j)*LS,FQW_TILE(l,N)*LS,FRACA(l,N),RESFS(l,N), &
               RESFT(l,N)
      print 24,l,FTL_1(i,j)*CP,FTL_TILE(l,N)*CP,FLAND(l),TILE_FRAC(l,N), &
               VFRAC_TILE(l,N)
      print 25,l,Z0H_TILE(l,N),Z0M_EFF_TILE(l,N)
      print 26,l,TSTAR_LAND(i,j),TSTAR_TILE(l,N),QSTAR_TILE(l,N), &
               CANHC_TILE(l,N)
      print 27,l,CD_TILE(l,N),CH_TILE(L,N),U_S_STD_TILE(l,N),V_S_TILE(l,N)
      print 28,l,TL_1(i,j),TS1(l),TSTAR(I,J),TSTAR_LAND(i,j),T_SOIL(l,1:4)
21    format(1x,'sfRADNET',i4,'z',f8.5,1x,f8.5,1x,f7.1,3f6.0)
22    format(1x,'sfRHOKH',i4,'z',f8.5,1x,f8.1,1x,f5.1,1x,f5.1,1x,f6.5)
23    format(1x,'sfFQW_1',i4,'z',f8.2,1x,f8.2,1x,f5.3,1x,f4.0,1x,f4.0)
24    format(1x,'sfFTL_1',i4,'z',f8.3,1x,f8.3,1x,f5.3,1x,f6.3,1x,f6.2)
25    format(1x,'sfZ0M',i4,'z',f8.5,1x,f8.5,1x,f5.2)
26    format(1x,'sfTSTAR',i4,'z',f6.1,1x,f6.1,1x,f7.5,1x,f7.0)
27    format(1x,'sfCDCHt',i4,'z',f7.4,1x,f7.4,1x,f7.4,1x,f7.4)
28    format(1x,'sfTSOIL',i4,'z',10f6.1)
!
!
!     l=1000 has i=13,j=29
      print *,'l=2132,i=28,j=57'
      l=1000
      i=13
      j=29
      l=2132
      i=28
      j=57
      print 31,l,RHOKPM(l,1:5),RHOKPM(l,8),RHOKPM_POT(l,1:5),RHOKPM_POT(l,8), &
               RADNET_TILE(l,1:5),RADNET_TILE(l,8),SW_DOWN(i,j),SW_TILE(l,1), &
               SW_TILE(l,8),LW_DOWN(i,j)
      print 32,l,RHOKH_1(l,1:5),RHOKH_1(l,8), &
               SNOW_TILE(l,1:5),SNOW_TILE(l,8),TL_1(i,j),TS1(l),QS1(i,j)
      print 33,l,FQW_1(i,j)*LC,FQW_TILE(l,1:5)*LC,FQW_TILE(l,8)*LC, &
               FRACA(l,1:5),FRACA(l,8),RESFS(l,1:5),RESFS(l,8), &
               RESFT(l,1:5),RESFT(l,8)
      print 34,l,FTL_1(i,j)*CP,FTL_TILE(l,1:5)*CP,FTL_TILE(l,8)*CP, &
               FLAND(l),TILE_FRAC(l,1:5),TILE_FRAC(l,8), &
               VFRAC_TILE(l,1:5)
      print 35,l,Z0H_TILE(l,1:5),Z0H_TILE(l,8), &
               Z0M_EFF_TILE(l,1:5),Z0M_EFF_TILE(l,8)
      print 36,l,TSTAR_LAND(i,j),TSTAR_TILE(l,1:5),TSTAR_TILE(l,8), &
               QSTAR_TILE(l,1:5),QSTAR_TILE(l,8), &
               CANHC_TILE(l,1:5),CANHC_TILE(l,8)
      print 37,l,CD_TILE(l,1:5),CD_TILE(l,8),CH_TILE(L,1:5),CH_TILE(l,8), &
        U_S_STD_TILE(l,1:5),U_S_STD_TILE(l,8),V_S_TILE(l,1:5),V_S_TILE(l,8)
      print 38,l,TL_1(i,j),TS1(l),TSTAR(I,J),TSTAR_LAND(i,j),T_SOIL(l,1:4), &
               LAI_FT(l,1:5),VFRAC_TILE(l,1:5)
31    format(1x,'sfRADNET',i4,'z',6f6.5,1x,6f6.5,1x,6f5.0,1x,4f6.0)
32    format(1x,'sfRHOKH',i4,'z',6f6.5,1x,6f6.0,1x,f5.1,1x,f5.1,1x,f6.5)
33    format(1x,'sfFQW_1',i4,'z',f9.3,1x,6f7.2,1x,6f5.3,1x,6f4.0,1x,6f4.0)
34    format(1x,'sfFTL_1',i4,'z',f8.3,1x,6f6.1,1x,f5.3,1x,6f5.3,1x,6f4.2)
35    format(1x,'sfZ0M',i4,'z',6f7.5,1x,6f7.5,1x,6f4.2)
36    format(1x,'sfTSTAR',i4,'z',f6.1,1x,6f5.0,1x,6f7.5,1x,6f7.0)
37    format(1x,'sfCDCHt',i4,'z',6f7.4,1x,6f7.4,1x,6f7.4,1x,6f7.4)
38    format(1x,'sfTSOIL,LAI',i4,'z',10f6.1,2x,5f4.2,1x,5f4.2)
!     l=1500 has i=94,j=44
      print *,'l=2080,i=35,j=56'
      l=1500
      i=94 
      j=44
      l=2080
      i=35
      j=56
      print 31,l,RHOKPM(l,1:5),RHOKPM(l,8),RHOKPM_POT(l,1:5),RHOKPM_POT(l,8), &
               RADNET_TILE(l,1:5),RADNET_TILE(l,8),SW_DOWN(i,j),SW_TILE(l,1), &
               SW_TILE(l,8),LW_DOWN(i,j)
      print 32,l,RHOKH_1(l,1:5),RHOKH_1(l,8), &
               SNOW_TILE(l,1:5),SNOW_TILE(l,8),TL_1(i,j),TS1(l),QS1(i,j)
      print 33,l,FQW_1(i,j)*LC,FQW_TILE(l,1:5)*LC,FQW_TILE(l,8)*LC, &
               FRACA(l,1:5),FRACA(l,8),RESFS(l,1:5),RESFS(l,8), &
               RESFT(l,1:5),RESFT(l,8)
      print 34,l,FTL_1(i,j)*CP,FTL_TILE(l,1:5)*CP,FTL_TILE(l,8)*CP, &
               FLAND(l),TILE_FRAC(l,1:5),TILE_FRAC(l,8), &
               VFRAC_TILE(l,1:5)
      print 35,l,Z0H_TILE(l,1:5),Z0H_TILE(l,8), &
               Z0M_EFF_TILE(l,1:5),Z0M_EFF_TILE(l,8)
      print 36,l,TSTAR_LAND(i,j),TSTAR_TILE(l,1:5),TSTAR_TILE(l,8), &
               QSTAR_TILE(l,1:5),QSTAR_TILE(l,8), &
               CANHC_TILE(l,1:5),CANHC_TILE(l,8)
      print 37,l,CD_TILE(l,1:5),CD_TILE(l,8),CH_TILE(L,1:5),CH_TILE(l,8), &
        U_S_STD_TILE(l,1:5),U_S_STD_TILE(l,8),V_S_TILE(l,1:5),V_S_TILE(l,8)
      print 38,l,TL_1(i,j),TS1(l),TSTAR(I,J),TSTAR_LAND(i,j),T_SOIL(l,1:4), &
               LAI_FT(l,1:5),VFRAC_TILE(l,1:5)
!      print 11,l,RHOKPM(l,:),RHOKPM_POT(l,:),RADNET_TILE(l,:)
!      print 12,l,RHOKH_1(l,:),SNOW_TILE(l,:),TL_1(i,j),TS1(l),QS1(i,j)
!      print 13,l,FQW_1(i,j),FQW_TILE(l,:),FRACA(l,:),RESFS(l,:) &
!                ,RESFT(l,:)
!      print 14,l,FTL_1(i,j),FTL_TILE(l,:),FLAND(l),TILE_FRAC(l,:) &
!                 ,VFRAC_TILE(l,:)
!      print 15,l,Z0H_TILE(l,:),Z0M_EFF_TILE(l,:)
!      print 16,l,TSTAR_LAND(i,j),TSTAR_TILE(l,:),QSTAR_TILE(l,:), &
!               CANHC_TILE(l,:)
11    format(1x,'sfRADNET',i4,'z',9f5.4,1x,9f5.4,1x,9f5.0)
12    format(1x,'sfRHOKH',i4,'z',9f5.4,1x,9f6.0,1x,f5.1,1x,f5.1,1x,f6.5)
13    format(1x,'sfFQW_1',i4,'z',f9.2,1x,9f7.2,1x,9f5.3,1x,9f4.0,1x,9f4.0)
14    format(1x,'sfFTL_1',i4,'z',f8.3,1x,9f5.0,1x,f5.3,1x,9f5.3,1x,9f4.2)
15    format(1x,'sfZ0M',i4,'z',9f7.5,1x,9f7.5,1x,9f4.2)
16    format(1x,'sfTSTAR',i4,'z',f6.1,1x,9f5.0,1x,9f6.5,1x,9f5.0)
      ENDIF

! Adjust ASHTF for sens. heat flux to ground beneath coniferous canopy
      IF( .NOT.  l_cable ) THEN
       IF (NTILES >  1 .AND. CAN_MODEL == 4) THEN
        N = 2
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          RH_CAN(L) = 43. / (SQRT(CD_TILE(L,N))*VSHR_LAND(I,J))
          ASHTF_TILE(L,N) = ASHTF_TILE(L,N) +                           &
     &                      VFRAC_TILE(L,N)*RHOSTAR(I,J)*CP/RH_CAN(L)
        ENDDO
       ENDIF
      ENDIF    ! (  .NOT. l_cable ) THEN

!-----------------------------------------------------------------------
!!  4.4   Calculate the standard deviations of layer 1 turbulent
!!        fluctuations of temperature and humidity using approximate
!!        formulae from first order closure.
!-----------------------------------------------------------------------

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        Q1_SD(I,J) = 0.
        T1_SD(I,J) = 0.
       ENDDO
      ENDDO

! Sea and sea-ice
! DEPENDS ON: stdev1_sea
      CALL STDEV1_SEA (                                                 &
     & ROW_LENGTH,ROWS,OFF_X,OFF_Y,FLANDG,                              &
     & BQ_1,BT_1,FQW_1,FTL_1,ICE_FRACT,RHOKM_SSI,RHOSTAR,VSHR_SSI,      &
     & Z0MSEA,Z0_ICE,Z1_TQ_SEA,                                         &
     & Q1_SD,T1_SD,LTIMER                                               &
     & )

! Land tiles
      DO N=1,NTILES
! DEPENDS ON: stdev1_land
        CALL STDEV1_LAND (                                              &
     &   ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS(N),                          &
     &   LAND_INDEX,TILE_INDEX(1,N),FLAND,                              &
     &   BQ_1,BT_1,FQW_TILE(1,N),FTL_TILE(1,N),RHOKM_1_TILE(1,N),       &
     &   RHOSTAR,VSHR_LAND,Z0M_TILE(1,N),Z1_TQ,TILE_FRAC(1,N),          &
     &   Q1_SD,T1_SD,LTIMER                                             &
     &   )
      ENDDO

!-----------------------------------------------------------------------
!! Calculate scaling parameters required for new boundary layer scheme
!-----------------------------------------------------------------------

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
!       ! Recall that CD is the dimensionless 'surface exchange
!       ! coefficient', see (1.1.23) in the documentation,
!       ! ie. lower case c_D.

! COR_UST == 1  gives the correct friction velocity
! which contains the square of VSR_LAND(I,J)
! otherwise the old formula is used for
! recovering the HADGEM1 configuration

          IF (COR_UST == 1) THEN
           U_S(I,J) = SQRT(                                             &
     &       FLANDG(I,J)*CD_LAND(I,J)*VSHR_LAND(I,J)*VSHR_LAND(I,J)     &
     &       +(1.-FLANDG(I,J))*CD_SSI(I,J)*VSHR_SSI(I,J)*VSHR_SSI(I,J)  &
     &                 )
          ELSE
           U_S(I,J) = SQRT(FLANDG(I,J)*CD_LAND(I,J)*VSHR_LAND(I,J)      &
     &    +(1.-FLANDG(I,J))*CD_SSI(I,J)*VSHR_SSI(I,J))
          ENDIF

        FB_SURF(I,J) = G * ( BT_1(I,J)*FTL_1(I,J) +                     &
     &                     BQ_1(I,J)*FQW_1(I,J) ) / RHOSTAR(I,J)
! Obukhov length
        IF (BL_diag%L_oblen) THEN
          BL_diag%oblen(I,J)= -(U_S(I,J)*U_S(I,J)*U_S(I,J))             &
     &    /(VKMAN*FB_SURF(I,J))
        ENDIF
! Ustar
        IF (BL_diag%L_ustar) THEN
          BL_diag%ustar(I,J)=U_S(I,J)
        ENDIF
! Surface buoyancy flux
        IF (BL_diag%L_wbsurf) THEN
          BL_diag%wbsurf(I,J)=FB_SURF(I,J)
        ENDIF

       ENDDO
      ENDDO



!-----------------------------------------------------------------------
!!  4.6 For sea points, calculate the wind mixing energy flux and the
!!      sea-surface roughness length on the P-grid, using time-level n
!!      quantities.
!-----------------------------------------------------------------------

      DO J=1,ROWS
       DO I=1,ROW_LENGTH

        IF (SFME) FME(I,J) = 0.0

!        if (i==6 .and. j==1) then
!          print *, 'dhb17: FLANDG(I,J)= ',FLANDG(I,J)
!        endif

        IF (FLANDG(I,J) <  1.0) THEN

!        if (i==6 .and. j==1) then
!          print *, 'dhb18: RHOKM_SSI,VSHR_SSI= ',RHOKM_SSI(I,J),VSHR_SSI(I,J)
!          print *, 'dhb18: ICE_FRACT= ',ICE_FRACT(I,J) 
!        endif

          TAU = RHOKM_SSI(I,J) * VSHR_SSI(I,J)             ! P243.130

!        if (i==6 .and. j==1) then
!          print *, 'dhb18: TAU= ',TAU
!        endif

          IF (ICE_FRACT(I,J)  >   0.0)                                  &
     &      TAU = RHOSTAR(I,J) * CD_SEA(I,J)                            &
     &        * VSHR_SSI(I,J) * VSHR_SSI(I,J)

!        if (i==6 .and. j==1) then
!          print *, 'dhb19: RHOSTAR,CD_SEA= ',RHOSTAR(I,J), CD_SEA(I,J)
!          print *, 'dhb19: VSHR_SSI,VSHR_SSI= ',VSHR_SSI(I,J), VSHR_SSI(I,J)
!          print *, 'dhb19: TAU,RHOSEA= ',TAU,RHOSEA
!        endif

          IF (SFME)                                                     &
     &      FME(I,J) = (1.0-ICE_FRACT(I,J)) * TAU * SQRT(TAU/RHOSEA)
!                                                            ! P243.96
! Limit Z0MSEA to 0.154m for TAU very small

!       YMA: set the new Charnock expression
!        if (i==6 .and. j==1) then
!          print *, 'dhb20: RT4CK,VSHR_SSI,UMD4CK= ',RT4CK,VSHR_SSI(I,J),UMD4CK
!        endif
         XCK4Z0M=RT4CK*(VSHR_SSI(I,J)-UMD4CK)
         EXPOSI=EXP(XCK4Z0M)
         EXNEGA=1./EXPOSI   !EXP(-XCK4Z0M)
         CK4Z0M=((EXPOSI-EXNEGA)/(EXPOSI+EXNEGA))*0.5 + 0.5
!        if (i==6 .and. j==1) then
!          print *, 'dhb21: CK4Z0MH,CK4Z0ML= ',CK4Z0MH,CK4Z0ML
!          print *, 'dhb21: BF: CK4Z0M= ',CK4Z0M
!        endif
         CK4Z0M=(CK4Z0MH-CK4Z0ML)*CK4Z0M + CK4Z0ML
!        if (i==6 .and. j==1) then
!          print *, 'dhb21: AF: CK4Z0M= ',CK4Z0M
!        endif

!          Z0MSEA(I,j) = 1.54E-6 / (SQRT(TAU/RHOSTAR(I,J)) + 1.0E-5)     &
!     &               +  (CHARNOCK/G) * (TAU / RHOSTAR(I,J))

!        if (i==6 .and. j==1) then
!          print *, 'dhb22: TAU,RHOSTAR= ',TAU,RHOSTAR(I,J)
!          print *, 'dhb22: CK4Z0M,G= ', CK4Z0M,G
!        endif

          Z0MSEA(I,j) = 1.54E-6 / (SQRT(TAU/RHOSTAR(I,J)) + 1.0E-5)     &
     &               +  (CK4Z0M/G) * (TAU / RHOSTAR(I,J))

!        if (i==6 .and. j==1) then
!          print *, 'dhb23: BF: Z0HSEA,Z0MSEA= ',Z0HSEA, Z0MSEA(I,J)
!        endif

          Z0MSEA(I,J) = MAX ( Z0HSEA , Z0MSEA(I,J) )

!        if (i==6 .and. j==1) then
!          print *, 'dhb24: AF: Z0MSEA= ',Z0MSEA(I,J)
!        endif

!        if((i==6).and.(j==1)) then
!	print*,"WRONGYMA: i,j,Z0MSEA(I,J) =", i,j,Z0MSEA(I,J)
!	print*,"VSHR_SSI(I,J)=", VSHR_SSI(I,J)
!        print*,"RT4CK, UMD4CK, XCK4Z0M, CK4Z0M, CK4Z0ML, CK4Z0M=", &
!     &          RT4CK, UMD4CK, XCK4Z0M, CK4Z0M, CK4Z0ML, CK4Z0M
!        endif
!                                       ... P243.B6 (Charnock formula)
!                    TAU/RHOSTAR is "mod VS squared", see eqn P243.131
        ENDIF

       ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Calculate effective roughness lengths, orographic blending heights
! and gridbox-average Richardson numbers.
!-----------------------------------------------------------------------

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        RIB(I,J) = 0.
        Z0M_EFF(I,J) = 1.
        Z0H_EFF(I,J) = 0.

! Sea and sea-ice (leads ignored)
        IF ( .NOT.LAND_MASK(I,J) )H_BLEND_OROG(I,J) = H_BLEND_MIN
        IF (FLANDG(I,J) <  1.0) THEN
          RIB(I,J) = RIB_SEA(I,J)
          Z0M_EFF(I,J) = Z0MSEA(I,J)
          Z0M_GB(I,J)  = Z0MSEA(I,J)
          Z0H_EFF(I,J) = Z0H_SEA(I,J)
          IF ( ICE_FRACT(I,J)  >   0. ) THEN
            RIB(I,J) = RIB_ICE(I,J)
            Z0M_EFF(I,J) = Z0_ICE(I,J)
            Z0M_GB(I,J)  = Z0_ICE(I,J)
            Z0H_EFF(I,J) = Z0_ICE(I,J)
          ENDIF
          Z0MSSI(I,J)  = Z0M_EFF(I,J)
          Z0HSSI(I,J)  = Z0H_EFF(I,J)
        ENDIF

       ENDDO
      ENDDO

      IF ( L_spec_z0 ) THEN
! Check for prescribed surface roughness lengths specified in SCM
! NAMELIST.  If specified in the &INPROF then they will be used
! instead of Model calculated values
        DO J=1,ROWS
        DO I=1,ROW_LENGTH
          IF ( Z0M_SCM(i,j) > 0.0 ) THEN
            ! Set z0m from SCM namelist for sea points
            Z0MSEA(I,J)  = Z0M_SCM(I,J)
            Z0M_EFF(I,J) = Z0M_SCM(I,J)
            Z0M_GB(I,J)  = Z0M_SCM(I,J)
            Z0MSSI(I,J)  = Z0M_SCM(I,J)
          ENDIF
          IF ( Z0H_SCM(i,j) > 0.0 ) THEN
            ! Set z0h from SCM namelist for sea points
            Z0H_SEA(I,J) = Z0H_SCM(I,J)
            Z0H_EFF(I,J) = Z0H_SCM(I,J)
            Z0HSSI(I,J)  = Z0H_SCM(I,J)
          ENDIF
        ENDDO
        ENDDO
      ENDIF

      DO L = 1,LAND_PTS
        Z0M_GB_LAND(L) = 0.
        Z0H_GB_LAND(L) = 0.
        FZ0(L) = 0.
      ENDDO

! Copy sea/sea-ice roughness lengths and Richardson numbers onto the
! the land point tiles.

! Weight so that gridbox mean can be calculated:
      DO L=1,LAND_PTS
!CDIR NODEP
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        IF(FLAND(L) <  1.0)THEN
          RIB(I,J) = (1.0-FLAND(L))*RIB(I,J)
          FZ0(L) = (1.0-FLAND(L)) / (LOG(LB/Z0MSSI(I,J))**2)
        ENDIF
      ENDDO

      DO N=1,NTILES
!CDIR NODEP
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          RIB(I,J) = RIB(I,J) + FLAND(L)*TILE_FRAC(L,N)*RIB_TILE(L,N)
          FZ0(L) = FZ0(L)                                               &
     &      + FLAND(L)*TILE_FRAC(L,N) / (LOG(LB/Z0M_TILE(L,N))**2)
        ENDDO
      ENDDO
      DO L = 1,LAND_PTS
        Z0M_GB_LAND(L) = LB * EXP( - SQRT(1./FZ0(L)) )

        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        Z0M_GB(I,J) = Z0M_GB_LAND(L)
      ENDDO
!
!     ! Calculate grid-box mean z0h
!     !  - needed for Z0H_EFF, calculated in SF_OROG_GB
!
      IF (SZ0HEFF) THEN

        DO L = 1,LAND_PTS
          FZ0(L) = 0.
        ENDDO

        DO L=1,LAND_PTS
!CDIR NODEP
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          IF(FLAND(L) <  1.0)THEN
            FZ0(L) = (1.0-FLAND(L)) / (LOG(LB/Z0HSSI(I,J))**2)
          ENDIF
        ENDDO
        DO N=1,NTILES
!CDIR NODEP
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          FZ0(L) = FZ0(L)                                               &
     &      + FLAND(L)*TILE_FRAC(L,N) / (LOG(LB/Z0H_TILE(L,N))**2)
        ENDDO
        ENDDO
        DO L = 1,LAND_PTS
          Z0H_GB_LAND(L) = LB * EXP( - SQRT(1./FZ0(L)) )
        ENDDO

      ENDIF  ! Calculate Z0H_EFF diagnostic

! DEPENDS ON: sf_orog_gb
      CALL SF_OROG_GB(                                                  &
     & ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX,                             &
     & LAND_MASK,FORMDRAG,FD_stab_dep,OROG_DRAG_PARAM,                  &
     & HO2R2_OROG,RIB,SIL_OROG_LAND,Z0M_GB_LAND,Z1_UV,                  &
     & H_BLEND_OROG,Z0M_EFF,SZ0HEFF,Z0H_GB_LAND,Z0H_EFF,LTIMER          &
     & )

!-----------------------------------------------------------------------
! If sea ice is present then set RECIP_L_MO_SEA to its value over 
! the ice, a long-standing choice for the screen diagnostics.
! Note that RECIP_L_MO_SEA is also used in BDY_EXPL2 to diagnose 
! shear-dominated boundary layer types.  To preserve bit 
! reproducibility when screen diagnostics are switched on, 
! this change has been moved outside the if-test on stash logicals
!-----------------------------------------------------------------------
      DO J=1,ROWS
        DO I=1,ROW_LENGTH
          IF (FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0. ) THEN
            RECIP_L_MO_SEA(I,J) = RECIP_L_MO_ICE(I,J)
          END IF
        END DO
      END DO
!-----------------------------------------------------------------------
!! Call SFL_INT to calculate CDR10M and CHR1P5M - interpolation coeffs
!! used to calculate screen temperature, humidity and 10m winds.
!-----------------------------------------------------------------------

      IF (SU10 .OR. SV10 .OR. SQ1P5 .OR. ST1P5) THEN

! Sea and sea-ice (leads ignored)
        DO J=1,ROWS
        DO I=1,ROW_LENGTH
          CDR10M(I,J) =0.
          IF (FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) <= 0. )           &
     &      Z0M_SEA(I,J) = Z0MSEA(I,J)
          IF (FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0. ) THEN
            CD_SEA(I,J) = CD_ICE(I,J)
            CH_SEA(I,J) = CH_ICE(I,J)
            Z0H_SEA(I,J) = Z0_ICE(I,J)
            Z0M_SEA(I,J) = Z0_ICE(I,J)
            V_S_SEA(I,J) = V_S_ICE(I,J)
          END IF
        END DO
        END DO

! DEPENDS ON: sfl_int_sea
        CALL SFL_INT_SEA (                                              &
     &   ROW_LENGTH,ROWS,OFF_X,OFF_Y,FLANDG,                            &
     &   VSHR_SSI,CD_SEA,CH_SEA,Z0M_SEA,Z0H_SEA,I_SCRN_T_DIAG,          &
     &   RECIP_L_MO_SEA,V_S_SEA,                                        &
     &   ICE_FRACT,NSICE,SICE_INDEX,Z1_UV,Z1_TQ_SEA,DB_ICE,             &
     &   SU10,SV10,ST1P5,SQ1P5,                                         &
     &   CDR10M,CHR1P5M_SICE,LTIMER                                     &
     &   )

! Land tiles
        DO N=1,NTILES
! DEPENDS ON: sfl_int_land
          CALL SFL_INT_LAND (                                           &
     &     ROW_LENGTH,ROWS,OFF_X,OFF_Y,LAND_PTS,TILE_PTS(N),            &
     &     TILE_INDEX(1,N),LAND_INDEX,FLANDG,                           &
     &     VSHR_LAND,CD_STD(1,N),CD_TILE(1,N),CH_TILE(1,N),             &
     &     TILE_FRAC(1,N),                                              &
     &     Z0M_EFF_TILE(1,N),Z0M_TILE(1,N),Z0H_TILE(1,N),               &
     &     I_SCRN_T_DIAG,RECIP_L_MO_TILE(1,N),                          &
     &     V_S_TILE(1,N),V_S_STD(1,N),                                  &
     &     Z1_UV,Z1_TQ,DB_TILE(1,N),                                    &
     &     SU10,SV10,ST1P5,SQ1P5,                                       &
     &     CDR10M,CHR1P5M(1,N),LTIMER                                   &
     &     )
        END DO

      ENDIF  ! IF (SU10 .OR. SV10 .OR. SQ1P5 .OR. ST1P5)

      IF(ntest>0) THEN

      print *,'CD10M1',CDR10M(4,2),CHR1P5M(100,9),RECIP_L_MO_TILE(100,9)
      print *,'CDM2',CDR10M(28,57),CHR1P5M(2132,1),RECIP_L_MO_TILE(2132,1)
      print *,'CDM3',CDR10M(35,56),CHR1P5M(2080,1),RECIP_L_MO_TILE(2080,1)
      ENDIF

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFEXCH  ',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_EXCH
