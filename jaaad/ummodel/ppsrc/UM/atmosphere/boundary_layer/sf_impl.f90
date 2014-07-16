
!jhan:there are a lot of sections here which i suspect can be deleted, or 
!jhan:attributed to "ifdef cable"

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE SF_IMPL-----------------------------------------------
!!!
!!!  Purpose: Calculate implicit correction to surface fluxes of heat,
!!!           moisture and momentum. Also calculates screen level
!!!           temperature and humidity as well as 10 m winds.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!  5.4   28/08/02  Change call to SF_MELT to enable canopy snow
!!!                  for Needleleaf trees.  R. Essery
!!!  5.4   28/08/02 Bug fix for sea ice.  M. Best
!!!  5.5   07/02/03 Added ice catagories  J.Ridley
!!!  6.1  01/09/04  Pass potential evaporation related variables.
!                                                          Nic Gedney
!!!  6.2  21/03/05  Pass through implicit scheme weight and
!!!                 remove hardwiring.
!!!                                          M. Diamantakis
!   6.2  06/01/06 Make a LW correction for each seaice catagory
!                                                          J.Ridley
! 6.2      21/02/06    Switch for mixing ratios   A.P.Lock
!!!  6.2  02/02/06  Passes L_flux_bc through argument list to allow
!!!                 settings for prescribed surface flux forcing
!!!                                                        R.Wong
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version ?, dated ?.
!!!
!!!  System component covered: P24.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE SF_IMPL (                                              &

! IN values defining field dimensions and subset to be processed :
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,LAND_PTS,                     &

! IN soil/vegetation/land surface data :
     & LAND_INDEX,LAND_MASK,NICE,                                       &
     & NTILES,TILE_INDEX,TILE_PTS,SM_LEVELS,                            &
     & CANHC_TILE,CANOPY,FLAKE,SMC,                                     &
     & TILE_FRAC,WT_EXT_TILE,                                           &
     & FLAND,FLANDG,                                                    &

! IN sea/sea-ice data :
     & DI,ICE_FRACT,DI_NCAT,ICE_FRACT_NCAT,U_0,V_0,                     &

! IN everything not covered so far :
     & PSTAR,LW_DOWN,SURF_RADFLUX,SW_TILE,TIMESTEP,                         &
     & T_SOIL,QW_1,TL_1,U_1,V_1,RHOKM_U_1,RHOKM_V_1,GAMMA,              &
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,                             &
     & DTRDZ_CHARNEY_GRID_1,DU_1,DV_1,                                  &
     & FQW_TILE,EPOT_TILE,FQW_ICE,FTL_ICE,                              &
     & FRACA,RESFS,RESFT,RHOKH,RHOKH_TILE,RHOKH_SICE,                   &
     & RHOKPM,RHOKPM_POT,RHOKPM_SICE,Z1,                                &
     & Z0HSSI,Z0MSSI,Z0H_TILE,Z0M_TILE,CDR10M_U,CDR10M_V,               &
     & CHR1P5M,CHR1P5M_SICE,CT_CTQ_1,DQW_1,DTL_1,CQ_CM_U_1,CQ_CM_V_1,   &
     & L_NEG_TSTAR,                                                     &
     & FLANDG_U,FLANDG_V,ANTHROP_HEAT,L_SICE_HEATFLUX,                  &

!    EAK
!    IN
     &  l_cable                                                         &
     &, surf_down_sw,alb_tile,cos_zenith_angle,l_tile_pts               &
!     &, surf_down_sw,alb_tile,cos_zenith_angle               &
     &, lat,long,day,time_sec                                           &
     &, ls_rain, ls_snow, conv_rain, conv_snow                          &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS,clapp                              &
     &, SATHH,SATCON,HCAP,HCON,Z1_UV                                    &
     &, smvccl, smvcwt, smvcst, sthf, sthu                              &
     &, SMCL                                                            &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &, FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB                         &
     &, TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB                     &
     &, USTAR_CAB,SURF_HTF_CAB                                          &
!
     &, TOT_ALB                                                         &  
     &, SNAGE_TILE,RTSOIL_TILE                                          &
     &, GFLUX_TILE,SGFLUX_TILE                                          &
     &, GC,GS,CANOPY_GB                                                 & ! Added GC, pfv 25oct13a
! Lestevens March 2010
     &, DIM_CS1, DIM_CS2    &
     &, NPP, NPP_FT         &
     &, GPP, GPP_FT         &
     &, RESP_S, RESP_S_TOT  &
     &, RESP_S_TILE         &  !kdcorbin, 10/10
     &, RESP_P, RESP_P_FT   &
     &, G_LEAF              &  !kdcorbin, 10/10
     &, TRANSP_TILE,        &  !Lestevens 3Nov11
! Lestevens Sept2012: CasaCNP variables 
     &  CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,                               &
     &  SOIL_ORDER,GLAI,PHENPHASE,                                      &
     ! Lestevens 23apr13
     &  NPP_FT_ACC,RESP_W_FT_ACC,                                        &

! IN STASH flags :-
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,                            &

! INOUT data :
     & TI,TI_GB,TSTAR,                                                  &
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,                       &
     & TSTAR_TILE,SNOW_TILE,                                            &
     & LE_TILE,RADNET_SICE,RADNET_TILE,                                 &
     & E_SEA,FQW_1,FTL_1,FTL_TILE,H_SEA,OLR,TAUX_1,TAUY_1,              &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &

! OUT Diagnostic not requiring STASH flags :
     & ECAN,EI_TILE,ESOIL_TILE,                                         &
     & SEA_ICE_HTF,SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,    &

! OUT diagnostic requiring STASH flags :
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,                        &
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE,U10M,V10M,                     &

! OUT data required elsewhere in UM system :
     & ECAN_TILE,EI,ESOIL,EXT,SNOWMELT,MELT_TILE,RHOKH_MIX,             &
     & ERROR,                                                           &

! LOGICAL LTIMER
     & lq_mix_bl,                                                       &
     & L_FLUX_BC,                                                       &
     & LTIMER                                                           &
     & )
 
      USE rad_switches_mod, ONLY: LRAD_EMIS_LAND_GEN,RAD_EMIS_LAND_GEN

      USE auscom_cpl_data_mod,                                          &
     &    Only : auscom_salinity, access_tfs, ocn_sss


      use cable_iday_mod
      IMPLICIT NONE

! EAK
      LOGICAL                                                           &
     & l_cable

!  Inputs :-
! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                   ! Local number of points on a row
     &, ROWS                                                            &
                   ! Local number of rows in a theta field
     &, N_ROWS                                                          &
                   ! Local number of rows in a v field
     &, OFF_X                                                           &
                   ! Size of small halo in i
     &, OFF_Y                                                           &
                   ! Size of small halo in j.
     &,LAND_PTS    ! IN No of land points

! (c) Soil/vegetation/land surface parameters (mostly constant).

      LOGICAL                                                           &
     & LAND_MASK(ROW_LENGTH,ROWS)  ! IN T if land, F elsewhere.

      INTEGER                                                           &
     & LAND_INDEX(LAND_PTS)        ! IN LAND_INDEX(I)=J => the Jth
!                                  !    point in ROW_LENGTH,ROWS is the
!                                  !    Ith land point.

      INTEGER                                                           &
     & SM_LEVELS                                                        &
                                   ! IN No. of soil moisture levels
     &,NTILES                                                           &
                                   ! IN No. of land tiles
     &,TILE_INDEX(LAND_PTS,NTILES)                                      &
                                   ! IN Index of tile points
     &,TILE_PTS(NTILES)                                                 &
                                   ! IN Number of tile points
     &,NICE                &       ! IN Number of sea ice catagories
! Lestevens March 2010
     &, DIM_CS1            &
     &, DIM_CS2

      REAL                                                              &
     & CANHC_TILE(LAND_PTS,NTILES)                                      &
                                   ! IN Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).
     &,CANOPY(LAND_PTS,NTILES)                                          &
                                   ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
     &,FLAKE(LAND_PTS,NTILES)                                           &
                                   ! IN Lake fraction.
     &,SMC(LAND_PTS)                                                    &
                                   ! IN Available soil moisture (kg/m2).
     &,TILE_FRAC(LAND_PTS,NTILES)                                       &
                                   ! IN Tile fractions including
!                                  ! snow cover in the ice tile.
     &,WT_EXT_TILE(LAND_PTS,SM_LEVELS,NTILES)                           &
!                                  ! IN Fraction of evapotranspiration
!                                  !    extracted from each soil layer
!                                  !    by each tile.
     &,FLAND(LAND_PTS)                                                  &
                                   ! IN Land fraction on land pts.
     &,FLANDG(ROW_LENGTH,ROWS)                                          &
!                                  ! IN Land fraction on all pts.
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
                                   ! IN Land fraction on U grid.
     &,FLANDG_V(ROW_LENGTH,N_ROWS) ! IN Land fraction on V grid.

! (d) Sea/sea-ice data.

      REAL                                                              &
     & DI(ROW_LENGTH,ROWS)                                              &
                                   ! IN "Equivalent thickness" of
!                                  !     sea-ice GBM agregate (m).
     &,DI_NCAT(ROW_LENGTH,ROWS,NICE)                                    &
                                     ! IN "Equivalent thickness" of
!                                  !     sea-ice catagories (m).
     &,ICE_FRACT(ROW_LENGTH,ROWS)                                       &
                                   ! IN Fraction of gridbox covered by
!                                  !     sea-ice (decimal fraction).
     &,ICE_FRACT_NCAT(ROW_LENGTH,ROWS,NICE)                             &
                                             ! IN Fraction of gridbox
!                                  !  covered by sea-ice on catagories.

     &,U_0(ROW_LENGTH,ROWS)                                             & 
                                   ! IN W'ly component of surface
!                                  !    current (m/s).
     &,V_0(ROW_LENGTH,N_ROWS)      ! IN S'ly component of surface
!                                  !    current (m/s).

! (f) Atmospheric + any other data not covered so far, incl control.

      REAL                                                              &
     & PSTAR(ROW_LENGTH,ROWS)                                           &
                                   ! IN Surface pressure (Pascals).
     &,LW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface downward LW radiation
!                                  !    (W/m2).
     &,SURF_RADFLUX(ROW_LENGTH,ROWS)                                        &
                                   ! IN Surface net SW and downward LW
!                                  !    radiation for sea-ice (W/sq m).
     &,SW_TILE(LAND_PTS,NTILES)                                         &
                                   ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
     &,TIMESTEP                                                         &
                                   ! IN Timestep (seconds).
     &,T_SOIL(LAND_PTS,SM_LEVELS)                                       &
                                   ! IN Soil temperatures (K).
     &,QW_1(ROW_LENGTH,ROWS)                                            &
                                   ! IN Total water content
     &,TL_1(ROW_LENGTH,ROWS)                                            &
                                   ! IN Ice/liquid water temperature
     &,U_1(1-OFF_X:ROW_LENGTH+OFF_X,                                    &
     &     1-OFF_Y:ROWS+OFF_Y)                                          & 
                                   ! IN W'ly wind component (m/s)
     &,V_1(1-OFF_X:ROW_LENGTH+OFF_X,                                    &
     &     1-OFF_Y:N_ROWS+OFF_Y)                                        & 
                                   ! IN S'ly wind component (m/s)
     &,RHOKM_U_1(ROW_LENGTH,ROWS)                                       &
                                   ! IN Exchange coefficients for
!                                  !    momentum (on U-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
     &,RHOKM_U_LAND(ROW_LENGTH,ROWS)                                    &
!                                  ! IN Exchange coefficients for
!                                  !    land momentum
!                                  !    (on U-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
     &,RHOKM_U_SSI(ROW_LENGTH,ROWS)                                     &
                                   ! IN Exchange coefficients for
!                                  !    mean sea momentum
!                                  !    (on U-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
     &,RHOKM_V_1(ROW_LENGTH,N_ROWS)                                     &
                                   ! IN Exchange coefficients for
!                                  !    momentum (on V-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
     &,RHOKM_V_LAND(ROW_LENGTH,N_ROWS)                                  &
!                                  ! IN Exchange coefficients for
!                                  !    land momentum
!                                  !    (on V-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
     &,RHOKM_V_SSI(ROW_LENGTH,N_ROWS)
!                                  ! IN Exchange coefficients for
!                                  !    mean sea momentum
!                                  !    (on V-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")

      REAL GAMMA                   ! IN implicit weight in level 1
      REAL                                                              &
     & ALPHA1(LAND_PTS,NTILES)                                          &
                                   ! IN Mean gradient of saturated
!                                  !    specific humidity with respect
!                                  !    to temperature between the
!                                  !    bottom model layer and tile
!                                  !    surfaces
     &,ALPHA1_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! IN ALPHA1 for sea-ice.
     &,ASHTF(ROW_LENGTH,ROWS)                                           &
                                   ! IN Coefficient to calculate surface
!                                  !    heat flux into soil or sea-ice.
     &,ASHTF_TILE(LAND_PTS,NTILES)                                      &
                                   !IN Coefficient to calculate
!                                  !    surface heat flux into land
!                                  !    tiles.
     &,DTRDZ_CHARNEY_GRID_1(ROW_LENGTH,ROWS)                            &
!                                  ! IN -g.dt/dp for model layers.
     &,FQW_TILE(LAND_PTS,NTILES)                                        &
                                   ! IN Surface FQW for land tiles
     &,EPOT_TILE(land_pts,ntiles)                                       &
                                   ! IN surface tile potential
!                                  !    evaporation
     &,FQW_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface FQW for sea-ice
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface FTL for sea-ice
     &,FRACA(LAND_PTS,NTILES)                                           &
                                   ! IN Fraction of surface moisture
!                                  !    flux with only aerodynamic
!                                  !    resistance for snow-free land
!                                  !    tiles.
     &,RESFS(LAND_PTS,NTILES)                                           &
                                   ! IN Combined soil, stomatal
!                                  !    and aerodynamic resistance
!                                  !    factor for fraction (1-FRACA) of
!                                  !    snow-free land tiles.
     &,RESFT(LAND_PTS,NTILES)                                           &
                                   ! IN Total resistance factor.
!                                  !    FRACA+(1-FRACA)*RESFS for
!                                  !    snow-free land, 1 for snow.
     &,RHOKH(ROW_LENGTH,ROWS)                                           &
                                   ! IN Grid-box surface exchange
!                                  !     coefficients
     &,RHOKH_TILE(LAND_PTS,NTILES)                                      &
                                   ! IN Surface exchange coefficients
!                                  !    for land tiles
     &,RHOKH_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! IN Surface exchange coefficients
!                                  !    for sea and sea-ice
     &,RHOKPM(LAND_PTS,NTILES)                                          &
                                   ! IN Land surface exchange coeff.
     &,RHOKPM_POT(LAND_PTS,NTILES)                                      &
                                   ! IN Land surface exchange coeff.
!                                    for potential evaporation.
     &,RHOKPM_SICE(ROW_LENGTH,ROWS)! IN Sea-ice surface exchange coeff.

       REAL                                                             &
     & Z1(ROW_LENGTH,ROWS)                                              &
                                   ! IN Height of lowest level (i.e.
!                                  !    height of middle of lowest
!                                  !    layer).
     &,Z0HSSI(ROW_LENGTH,ROWS)                                          &
     &,Z0MSSI(ROW_LENGTH,ROWS)                                          &
                                   ! IN Roughness lengths over sea (m)
     &,Z0H_TILE(LAND_PTS,NTILES)                                        &
                                   ! IN Tile roughness lengths for heat
!                                  !    and moisture (m).
     &,Z0M_TILE(LAND_PTS,NTILES)                                        &
                                   ! IN Tile roughness lengths for
!                                  !    momentum.
     &,CDR10M_U(ROW_LENGTH,ROWS)                                        & 
                                   ! IN Ratio of CD's reqd for
!                                  !    calculation of 10 m wind. On
!                                  !    U-grid; comments as per RHOKM.
     &,CDR10M_V(ROW_LENGTH,N_ROWS)                                      & 
                                   ! IN Ratio of CD's reqd for
!                                  !    calculation of 10 m wind. On
!                                  !    V-grid; comments as per RHOKM.
     &,CHR1P5M(LAND_PTS,NTILES)                                         &
                                   ! IN Ratio of coefffs for calculation
!                                  !    of 1.5m temp for land tiles.
     &,CHR1P5M_SICE(ROW_LENGTH,ROWS)                                    &
!                                  ! IN CHR1P5M for sea and sea-ice
!                                  !    (leads ignored).
     &,CT_CTQ_1(ROW_LENGTH,ROWS)                                        &
                                   ! IN Coefficient in T and q
!                                  !    tri-diagonal implicit matrix
     &,CQ_CM_U_1(ROW_LENGTH,ROWS)                                       &
                                   ! IN Coefficient in U tri-diagonal
!                                  !    implicit matrix
     &,CQ_CM_V_1(ROW_LENGTH,N_ROWS)                                     &
                                   ! IN Coefficient in V tri-diagonal
!                                  !    implicit matrix
     &,DQW_1(ROW_LENGTH,ROWS)                                           &
                                   ! IN Level 1 increment to q field
     &,DTL_1(ROW_LENGTH,ROWS)                                           &
                                   ! IN Level 1 increment to T field
     &,DU_1(1-OFF_X:ROW_LENGTH+OFF_X,                                   &
     &      1-OFF_Y:ROWS+OFF_Y)                                         &
                                   ! IN Level 1 increment to u wind
!                                  !    field
     &,DV_1(1-OFF_X:ROW_LENGTH+OFF_X,                                   &
     &      1-OFF_Y:N_ROWS+OFF_Y)                                       &
                                   ! IN Level 1 increment to v wind
!                                  !    field
     &,ANTHROP_HEAT(NTILES) 
                                   ! IN Additional heat source on tiles
                                   !    for anthropogenic urban heat
                                   !    source (W/m2)

      LOGICAL                                                           &
     & LTIMER                                                           &
                                   ! IN Logical switch for TIMER diags
     &,L_NEG_TSTAR                                                      &
                                   ! IN Switch for -ve TSTAR error check
     &,L_FLUX_BC                                                        &
                                   ! IN SCM logical for prescribed
                                   !    surface flux forcing
     &,L_SICE_HEATFLUX             ! IN T: semi-implicit update of TI

      LOGICAL                                                           &
     & lq_mix_bl              ! TRUE if mixing ratios used in
!                             ! boundary layer code

!  STASH flags :-

      LOGICAL                                                           &
     & SIMLT                                                            &
               ! IN Flag for SICE_MLT_HTF (q.v.)
     &,SMLT                                                             &
               ! IN Flag for SNOMLT_SURF_HTF (q.v.)
     &,SLH                                                              &
               ! IN Flag for LATENT_HEAT (q.v.)
     &,SQ1P5                                                            &
               ! IN Flag for Q1P5M (q.v.)
     &,ST1P5                                                            &
               ! IN Flag for T1P5M (q.v.)
     &,SU10                                                             &
               ! IN Flag for U10M (q.v.)
     &,SV10    ! IN Flag for V10M (q.v.)

!  In/outs :-

      REAL                                                              &
     & TI(ROW_LENGTH,ROWS,NICE)                                         &
                                   ! INOUT Sea-ice surface layer
!                                  !       temperature (K).
     &,TI_GB(ROW_LENGTH,ROWS)                                           &
                                   ! OUT GBM ice surface temperature (K)
     &,TSTAR(ROW_LENGTH,ROWS)                                           &
                                   ! OUT   GBM surface temperature (K).
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! OUT   Land mean sfc temperature (K)
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! IN    Open sea sfc temperature (K).
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! OUT   Sea-ice sfc temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)                                       &
                                   ! INOUT Sea mean sfc temperature (K).
     &,TSTAR_TILE(LAND_PTS,NTILES)                                      &
                                   ! INOUT Surface tile temperatures
     &,SNOW_TILE(LAND_PTS,NTILES)                                       &
                                   ! INOUT Lying snow on tiles (kg/m2)
     &,LE_TILE(LAND_PTS,NTILES)                                         &
                                   ! INOUT Surface latent heat flux for
!                                  !     land tiles
     &,RADNET_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! INOUT Surface net radiation on
!                                  !       sea-ice (W/m2)
     &,RADNET_TILE(LAND_PTS,NTILES)                                     &
                                   ! INOUT Surface net radiation on
!                                  !       land tiles (W/m2)
     &,E_SEA(ROW_LENGTH,ROWS)                                           &
                                   ! INOUT Evaporation from sea times
!                                  !       leads fraction. Zero over
!                                  !       land. (kg per square metre
!                                  !       per sec).
     &,FQW_1(ROW_LENGTH,ROWS)                                           &
                                   ! INOUT Moisture flux between layers
!                                  !       (kg per square metre per sec)
!                                  !       FQW(,1) is total water flux
!                                  !       from surface, 'E'.
     &,FTL_1(ROW_LENGTH,ROWS)                                           &
                                   ! INOUT FTL(,K) contains net
!                                  !       turbulent sensible heat flux
!                                  !       into layer K from below; so
!                                  !       FTL(,1) is the surface
!                                  !       sensible heat, H.(W/m2)
     &,FTL_TILE(LAND_PTS,NTILES)                                        &
                                   ! INOUT Surface FTL for land tiles
     &,H_SEA(ROW_LENGTH,ROWS)                                           &
                                   ! INOUT Surface sensible heat flux
!                                  !       over sea times leads fraction
!                                  !       (W/m2)
     &,OLR(ROW_LENGTH,ROWS)                                             &
                                   ! IN    TOA - surface upward LW on
!                                  !       last radiation timestep
!                                  ! OUT   Corrected TOA outward LW
     &,TAUX_1(ROW_LENGTH,ROWS)                                          & 
                                   ! OUT   W'ly component of surface
!                                  !       wind stress (N/sq m). (On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or, at
!                                  !       present, set to missing data
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       & 
                                   ! INOUT W'ly component of surface
!                                  !       wind stress over land
!                                  !       (N/sq m). (On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or, at
!                                  !       present, set to missing data
     &,TAUX_SSI(ROW_LENGTH,ROWS)                                        & 
                                   ! INOUT W'ly component of surface
!                                  !       wind stress over mean sea
!                                  !       (N/sq m). (On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or, at
!                                  !       present, set to missing data
     &,TAUY_1(ROW_LENGTH,N_ROWS)                                        & 
                                   ! OUT   S'ly component of surface
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     & 
                                   ! INOUT S'ly component of land sfc
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX
     &,TAUY_SSI(ROW_LENGTH,N_ROWS) ! INOUT S'ly compt of mean sea sfc
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX

!     EAK
      REAL                                                              &
     &  alb_tile(land_pts,ntiles,4)                                  &
     &, surf_down_sw(row_length,rows,4)                                 &
     &, cos_zenith_angle(row_length,rows)                               &
     &, lat(row_length, rows)                                           &
                                          ! Lat. of gridpoint chosen
     &, long(row_length, rows)                                          &
                                          ! Long. of gridpoint chosen
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)                                     &
     &, smvccl (land_pts)                                            &
                             ! soil/qrparm.soil.crit
     &, smvcwt (land_pts)                                            &
                             ! soil/qrparm.soil.wilt
     &, smvcst (land_pts)                                            &
                             ! soil/qrparm.soil.satn
     &, sthf(land_pts,sm_levels)                                    &
                                ! IN Frozen soil moisture content of
                                !     each layer as a fraction of
                                !     saturation.
     &, sthu(land_pts,sm_levels)                                    &
                                ! IN Unfrozen soil moisture content
                                !    of each layer as a fraction of
                                !    saturation.
     &, time_sec                                                        &
                                ! actual time of day in secs.
     &, SNOW_DEPTH3L(land_pts,NTILES,3)                          &
     &, SNOW_MASS3L(land_pts,NTILES,3)                           &
     &, SNOW_COND(land_pts,NTILES,3)                             &
     &, SNOW_TMP3L(land_pts,NTILES,3)                            &
     &, SNOW_RHO3L(land_pts,NTILES,3)                            &
     &, SNOW_RHO1L(land_pts,NTILES)                              &
     &, SNAGE_TILE(land_pts,NTILES)                              &
     &, SMCL_TILE(land_pts,NTILES,sm_levels)                     &
     &, STHU_TILE(land_pts,NTILES,sm_levels)                     &
     &, STHF_TILE(land_pts,NTILES,sm_levels)                     &
     &, TSOIL_TILE(land_pts,NTILES,sm_levels)                    &
     &, T_SURF_TILE(land_pts,NTILES)                             &
     &, RTSOIL_TILE(land_pts,NTILES)                             &
     &, GFLUX_TILE(land_pts,NTILES)                              &
     &, SGFLUX_TILE(land_pts,NTILES)                             &
     &, Z1_UV(row_length,rows)                                          &
                                     ! Height of lowest u,v level.
     &,U_S(ROW_LENGTH,ROWS)                                             & 
                                   ! OUT Surface friction velocity (m/s)
     &, HCONS(land_pts)
       Real                                                              &
     &  clapp(land_pts)                                              &
                             !  qrparm.soil.bwag ?
!                               Clapp-Hornberger exponent.
     &, satcon(land_pts)                                             &
                                 !  qrparm.soil.satcon
     &, sathh(land_pts)                                              &
                             !  soil water suction
     &, hcon (land_pts)                                              &
                             ! soil/qrparm.soil.hcond
     &, hcap (land_pts)                                              &
                             ! soil/qrparm.soil.hcap
     &, SMCL(land_pts, sm_levels)
                             !  qrclim.smc_pm(lev).(month)
!  EAK    diagnostic variables for CABLE output
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


      Real                                                           &
     & FTL_TILE_CAB(LAND_PTS,NTILES)                                 &
     &,FTL_CAB(LAND_PTS)                                             &
     &,LE_TILE_CAB(LAND_PTS,NTILES)                                  &
     &,LE_CAB(LAND_PTS)                                              &
     &,TSTAR_TILE_CAB(LAND_PTS,NTILES)                               &
     &,TSTAR_CAB(LAND_PTS)                                           &
     &,SMCL_CAB(LAND_PTS,SM_LEVELS)                                  &
     &,TSOIL_CAB(LAND_PTS,SM_LEVELS)                                 &
     &,USTAR_CAB(LAND_PTS)                                           &
!
     &,TOT_ALB(LAND_PTS,NTILES)                                      &
     &,SURF_HTF_CAB(LAND_PTS)                                        &
     &,CANOPY_GB(LAND_PTS)                                           &
     &,GC(LAND_PTS,NTILES)                                           & ! Added GC, pfv 25oct13a
     &,GS(LAND_PTS)                                                  &
! Lestevens 17 Feb 2010 - passing CO2 fluxes
     !kdcorbin, 11/10 - changed from NPFT
     &,NPP_FT(LAND_PTS,NTILES)                                       &
     &,NPP(LAND_PTS)                                                 &
     !kdcorbin, 11/10 - changed from NPFT
     &,GPP_FT(LAND_PTS,NTILES)                                       &
     &,GPP(LAND_PTS)                                                 &
     &,RESP_S(LAND_PTS,DIM_CS1)                                      &
     &,RESP_S_TOT(DIM_CS2)                                           &
     &,RESP_S_TILE(LAND_PTS,NTILES)                      & !kdcorbin, 10/10
     &,RESP_P(LAND_PTS)                                              &
     !kdcorbin, 11/10 - changed from NPFT
     &,RESP_P_FT(LAND_PTS,NTILES)                                    &
     &,G_LEAF(LAND_PTS,NTILES)         &
     &,TRANSP_TILE(LAND_PTS,NTILES)    &
! Lestevens Sept2012: CasaCNP variables 
     &, CPOOL_TILE(LAND_PTS,NTILES,10) &
     &, NPOOL_TILE(LAND_PTS,NTILES,10) &
     &, PPOOL_TILE(LAND_PTS,NTILES,12) &
     &, SOIL_ORDER(LAND_PTS)           &
     &, GLAI(LAND_PTS,NTILES)          &
     &, PHENPHASE(LAND_PTS,NTILES)     &
     ! Lestevens 23apr13
     &, NPP_FT_ACC(land_pts,NTILES)    &
     &, RESP_W_FT_ACC(land_pts,NTILES)

      Integer                                                           &
     &  day                                                             &
     &, total_nsteps                                                    &
                                ! Total number of steps in run
     &, SOIL_TYPE(row_length,rows)                                      &
     &, VEG_TYPE(row_length,rows)                                       &
     &, ISNOW_FLG3L(land_pts,NTILES)

      Logical                                                           &
     &  l_tile_pts(land_pts,ntiles)


!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL                                                              &
     & ECAN(ROW_LENGTH,ROWS)                                            &
                                   ! OUT Gridbox mean evaporation from
!                                  !     canopy/surface store (kg/m2/s).
!                                  !     Zero over sea.
     &,ESOIL_TILE(LAND_PTS,NTILES)                                      &
                                   ! OUT ESOIL for snow-free land tiles
     &,SEA_ICE_HTF(ROW_LENGTH,ROWS,NICE)                                &
                                        ! OUT Heat flux through sea-ice
!                                  !     (W/m2, positive downwards).
     &,SURF_HT_FLUX(ROW_LENGTH,ROWS)                                    &
!                                  ! OUT Net downward heat flux at
!                                  !     surface over land and sea-ice
!                                  !     fraction of gridbox (W/m2).
     &,SURF_HT_FLUX_LAND(ROW_LENGTH,ROWS)                               &
!                                  ! OUT Net downward heat flux at
!                                  !     surface over land
!                                  !     fraction of gridbox (W/m2).
     &,SURF_HT_FLUX_SICE(ROW_LENGTH,ROWS)
!                                  ! OUT Net downward heat flux at
!                                  !     surface over sea-ice
!                                  !     fraction of gridbox (W/m2).

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

      REAL                                                              &
     & SICE_MLT_HTF(ROW_LENGTH,ROWS,NICE)                               &
!                                  ! OUT Heat flux due to melting of
!                                  !     sea-ice (Watts per sq metre).
     &,SNOMLT_SURF_HTF(ROW_LENGTH,ROWS)                                 &
!                                  ! OUT Heat flux required for surface
!                                  !     melting of snow (W/m2).
     &,LATENT_HEAT(ROW_LENGTH,ROWS)                                     &
                                   ! OUT Surface latent heat flux, +ve
!                                  !     upwards (Watts per sq m).
     &,Q1P5M(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Q at 1.5 m (kg water / kg air).
     &,Q1P5M_TILE(LAND_PTS,NTILES)                                      &
                                   ! OUT Q1P5M over land tiles.
     &,T1P5M(ROW_LENGTH,ROWS)                                           &
                                   ! OUT T at 1.5 m (K).
     &,U10M(ROW_LENGTH,ROWS)                                            &
                                   ! OUT U at 10 m (m per s).
     &,T1P5M_TILE(LAND_PTS,NTILES)                                      &
                                   ! OUT T1P5M over land tiles.
     &,V10M(ROW_LENGTH,N_ROWS)     ! OUT V at 10 m (m per s).

!-2 Genuinely output, needed by other atmospheric routines :-

      REAL                                                              &
     & EI(ROW_LENGTH,ROWS)                                              &
                                   ! OUT Sublimation from lying snow or
!                                  !     sea-ice (kg/m2/s).
     &,EI_LAND(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Sublimation from lying snow
!                                  !     (kg/m2/s).
     &,EI_SICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Sublimation from sea-ice
!                                  !     (kg/m2/s).
     &,EI_TILE(LAND_PTS,NTILES)                                         &
                                   ! OUT EI for land tiles.
     &,ECAN_TILE(LAND_PTS,NTILES)                                       &
                                   ! OUT ECAN for snow-free land tiles
     &,ESOIL(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Surface evapotranspiration
!                                  !     from soil moisture store
!                                  !     (kg/m2/s).
     &,EXT(LAND_PTS,SM_LEVELS)                                          &
                                   ! OUT Extraction of water from each
!                                  !     soil layer (kg/m2/s).
     &,SNOWMELT(ROW_LENGTH,ROWS)                                        &
                                   ! OUT Snowmelt (kg/m2/s).
     &,MELT_TILE(LAND_PTS,NTILES)                                       &
                                   ! OUT Snowmelt on land tiles (kg/m2/s
     &,RHOKH_MIX(ROW_LENGTH,ROWS)  ! OUT Exchange coeffs for moisture.

      INTEGER                                                           &
     & ERROR          ! OUT 0 - AOK;
!                     !     1 to 7  - bad grid definition detected;
!jhan:we dont need these
!      REAL                                                              &
!     & LAT_HT1,LAT_HT2,LAT_HT3

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL IM_SF_PT,SF_EVAP,SF_MELT,SCREEN_TQ,SICE_HTF
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
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
! CSIGMA start
      ! Stefan-Boltzmann constant (W/m**2/K**4).
      REAL, PARAMETER ::  SBCON=5.67E-8
! CSIGMA end
! CAOPTR Holds address pointers for atmosphere-to-ocean coupling fields
! required by SWAP_A2O and SWAP_O2A.
! 5.5 28/02/03 Add pointers for river outflow on atmos grid
! later remove runoff and ocenpts pointers  C.Bunton
! 6.2 24/02/06 Add pointers for 10m windspeed            J.Gunson
! 6.2 29/11/05 Add pointers for ice velocities on OCEAN grid  A.B.Keen
!
! 6.2 24/02/06 Add pointers for ocean DMS flux to atmosphere    J.Gunson
! Pointers needed by SWAP_A2O and SWAP_O2A

      INTEGER :: JO_SST    ! Sea-surface temperature on ocean grid
      INTEGER :: JO_UCURR  ! Surface zonal current on ocean grid
      INTEGER :: JA_AICE   ! Fractional ice conc. on atmos grid
      INTEGER :: JA_AICEN  ! Category frac ice conc. on atmos grid
      INTEGER :: JO_AICEN  ! Category frac ice conc. on ocean grid
      COMMON /AOPTR/JO_SST,JO_UCURR,JA_AICE,JA_AICEN,JO_AICEN

      ! Pointers needed by SWAP_A2O

      INTEGER:: JA_TAUX       ! Surface x-windstress on atmos grid
      INTEGER:: JO_TAUX       ! Surface x-windstress on ocean grid
      INTEGER:: JA_TAUY       ! Surface y-windstress on atmos grid
      INTEGER:: JO_TAUY       ! Surface y-windstress on ocean grid
      INTEGER:: JA_WINDMIX    ! Windmixing power on atmos grid
      INTEGER:: JO_WINDMIX    ! Windmixing power on ocean grid
      INTEGER:: JA_W10        ! 10m windspeed
      INTEGER:: JO_W10
      INTEGER:: JA_DDDL1_D1   ! Dust dry dep lyr 1
      INTEGER:: JA_DDDL1_D2
      INTEGER:: JA_DDDL1_D3
      INTEGER:: JA_DDDL1_D4
      INTEGER:: JA_DDDL1_D5
      INTEGER:: JA_DDDL1_D6
      INTEGER:: JA_DDDL2_D1   ! Dust dry dep lyr 2
      INTEGER:: JA_DDDL2_D2
      INTEGER:: JA_DDDL2_D3
      INTEGER:: JA_DDDL2_D4
      INTEGER:: JA_DDDL2_D5
      INTEGER:: JA_DDDL2_D6
      INTEGER:: JA_DWDLS_D1   ! Dust wet dep lrg scl rain
      INTEGER:: JA_DWDLS_D2
      INTEGER:: JA_DWDLS_D3
      INTEGER:: JA_DWDLS_D4
      INTEGER:: JA_DWDLS_D5
      INTEGER:: JA_DWDLS_D6
      INTEGER:: JA_DWDCN_D1   ! Dust wet dep convctv rain
      INTEGER:: JA_DWDCN_D2
      INTEGER:: JA_DWDCN_D3
      INTEGER:: JA_DWDCN_D4
      INTEGER:: JA_DWDCN_D5
      INTEGER:: JA_DWDCN_D6
      INTEGER:: JO_DDEPD1     ! Dust deposition : ocean
      INTEGER:: JO_DDEPD2
      INTEGER:: JO_DDEPD3
      INTEGER:: JO_DDEPD4
      INTEGER:: JO_DDEPD5
      INTEGER:: JO_DDEPD6
      INTEGER:: JA_SOLAR      ! Net downward SW at surf on atmos grid
      INTEGER:: JA_BLUE       ! Net blueband SW at surf on atmos grid
      INTEGER:: JO_BLUE       ! Net blueband SW at surf on ocean grid
      INTEGER:: JA_EVAP       ! Net evaporation over sea on atmos grid
      INTEGER:: JA_LONGWAVE   ! Net downward LW at surf on atmos grid
      INTEGER:: JA_SENSIBLE  ! Sensible heat flux over sea on atmos grid
      INTEGER:: JO_HEATFLUX   ! Non penetrative heatflux on ocean grid
      INTEGER:: JA_LSSNOW     ! Large-scale snowfall rate on atmos grid
      INTEGER:: JA_CVSNOW     ! Convective snowfall rate on atmos grid
      INTEGER:: JA_LSRAIN     ! Large-scale rainfall rate on atmos grid
      INTEGER:: JA_CVRAIN     ! Convective rainfall rate on atmos grid
      INTEGER:: JO_PMINUSE    ! Precipitation-evaporation on ocean grid
      INTEGER:: JA_SLOWRUNOFF ! Slow (sub-surface) runoff on atmos grid
      INTEGER:: JA_FASTRUNOFF ! Fast (surface) runoff on atmos grid
      INTEGER:: JA_RIVEROUT   ! Total river outflow on atmos grid
      INTEGER:: JA_OCENTPTS   ! Ocean entry point index to atmos landpts
      INTEGER:: JO_RIVEROUT   ! Total river outflow on ocean grid
      INTEGER:: JA_co2        ! atmos level 1 co2 conc.
      INTEGER:: JO_co2
      INTEGER:: JA_co2flux    ! ocean co2 flux.
      INTEGER:: JO_co2flux
      INTEGER:: JA_dmsflux    ! ocean DMS flux.
      INTEGER:: JO_dmsflux
      INTEGER:: JO_SNOWFALL   ! Snowfall rate on ocean grid
      INTEGER:: JA_SUBLIM     ! Sublimation on atmos grid
      INTEGER:: JO_SUBLIM     ! Sublimation on ocean grid
      INTEGER:: JA_BOTMELT    ! Diffusive heat thro ice on atmos grid
      INTEGER:: JA_BOTMELTN   ! Diff heat thro ice (ncat)
      INTEGER:: JO_BOTMELT    ! Diffusive heat thro ice on ocean grid
      INTEGER:: JA_TOPMELT    ! Seaice top melting flux on atmos grid
      INTEGER:: JA_TOPMELTN   ! Seaice top melting flux (ncat)
      INTEGER:: JO_TOPMELT    ! Seaice top melting flux on ocean grid
      INTEGER:: JO_AICE       ! Sea ice concentration on ocean grid
      INTEGER:: JA_PRESS      ! Surface pressure on atmos grid
      INTEGER :: JC_U10
      INTEGER :: JC_V10
      INTEGER :: JC_CO2
 
      COMMON /A2OPTR/                                                   &
     &  JA_TAUX,JO_TAUX,JA_TAUY,JO_TAUY,JA_WINDMIX,JO_WINDMIX,          &
     &  JA_W10, JO_W10,                                                 &
     &  JA_DDDL1_D1,JA_DDDL1_D2,JA_DDDL1_D3,JA_DDDL1_D4,JA_DDDL1_D5,    &
     &  JA_DDDL1_D6,JA_DDDL2_D1,JA_DDDL2_D2,JA_DDDL2_D3,JA_DDDL2_D4,    &
     &  JA_DDDL2_D5,JA_DDDL2_D6,JA_DWDLS_D1,JA_DWDLS_D2,JA_DWDLS_D3,    &
     &  JA_DWDLS_D4,JA_DWDLS_D5,JA_DWDLS_D6,JA_DWDCN_D1,JA_DWDCN_D2,    &
     &  JA_DWDCN_D3,JA_DWDCN_D4,JA_DWDCN_D5,JA_DWDCN_D6,                &
     &  JO_DDEPD1,JO_DDEPD2,JO_DDEPD3,JO_DDEPD4,JO_DDEPD5,JO_DDEPD6,    &
     &  JA_SOLAR,JA_BLUE,JO_BLUE,JA_EVAP,JA_LONGWAVE,JA_SENSIBLE,       &
     &  JO_HEATFLUX,JA_LSSNOW,JA_CVSNOW,JA_LSRAIN,JA_CVRAIN,JO_PMINUSE, &
     &  JA_SLOWRUNOFF,JA_FASTRUNOFF,JA_OCENTPTS,JA_RIVEROUT,            &
     &  JO_RIVEROUT,JA_co2, JO_co2, JA_co2flux, JO_co2flux,             &
     &  JA_dmsflux, JO_dmsflux,                                         &
     &  JO_SNOWFALL,JA_SUBLIM,JO_SUBLIM,JA_BOTMELT,JO_BOTMELT,          &
     &  JA_TOPMELT,JO_TOPMELT,JA_BOTMELTN,JA_TOPMELTN,JO_AICE,JA_PRESS, &
     &  JC_U10, JC_V10, JC_CO2 

      INTEGER:: JO_TSTAR           ! Surface temperature on ocean grid
      INTEGER:: JA_TSTAR           ! Surface temperature on atmos grid
      INTEGER:: JA_UCURR           ! Surface zonal current on atmos grid
      INTEGER:: JO_VCURR           ! Surface merid current on ocean grid
      INTEGER:: JA_VCURR           ! Surface merid current on atmos grid
      INTEGER:: JO_UICE            ! Zonal seaice vel on ocean grid
      INTEGER:: JO_VICE            ! Merid seaice vel on ocean grid
      INTEGER:: JO_ICEDEPTH        ! Ice depth on ocean grid
      INTEGER:: JA_ICEDEPTH        ! Ice depth on atmos grid
      INTEGER:: JA_ICEDEPTHN       ! Ice depth on atmos grid (ncat)
      INTEGER:: JO_SNOWDEPTH       ! Snow depth on ocean grid
      INTEGER:: JA_SNOWDEPTH       ! Snow depth on atmos grid
      INTEGER:: JA_SNOWDEPTHN      ! Snow depth on atmos grid (ncat)

      COMMON /O2APTR/                                                   &
     &  JO_TSTAR,JA_TSTAR,JA_UCURR,JO_VCURR,JA_VCURR,                   &
     &  JO_UICE,JO_VICE,                                                &
     &  JO_ICEDEPTH,JA_ICEDEPTH,JO_SNOWDEPTH,JA_SNOWDEPTH,              &
     &  JA_ICEDEPTHN,JA_SNOWDEPTHN
! CAOPTR end
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

! Derived local parameters.

      REAL LS

      PARAMETER (                                                       &
     & LS=LF+LC                                                         &
                             ! Latent heat of sublimation.
     &  )

!-----------------------------------------------------------------------

!  Workspace :-

      REAL                                                              &
     & ELAKE_TILE(LAND_PTS,NTILES)                                      &
                                   ! Lake evaporation.
     &,QIM_1(ROW_LENGTH,ROWS)                                           &
                                   ! Implicit value of first model level
!                                  ! humidity
     &,TIM_1(ROW_LENGTH,ROWS)                                           &
                                   ! Implicit value of first model level
!                                  ! temperature
     &,TSTAR_RAD4(ROW_LENGTH,ROWS)                                      &
                                   ! Effective surface radiative
!                                  ! temperature for land and sea-ice
     &,TSTAR_TILE_OLD(LAND_PTS,NTILES)                                  &
!                                  ! Tile surface temperatures at
!                                  ! beginning of timestep.
     &,SICE_MELT(ROW_LENGTH,ROWS,NICE)                                  &
                                       !Melt at surface sea-ice catagory
     &,TSTAR_SIC(ROW_LENGTH,ROWS,NICE)                                  &
                                       !Ice catagory surface temperature
     &,SURF_HT_FLUX_SICE_NCAT(ROW_LENGTH,ROWS,NICE)
!                                  ! heat flux by ice catagory
      REAL TSTAR_SSI0(ROW_LENGTH,ROWS)
      REAL SICE_MELT0(ROW_LENGTH,ROWS)


!  Local scalars :-

      REAL                                                              &
     & LAT_HT     ! Latent heat of evaporation for snow-free land
!                 ! or sublimation for snow-covered land and ice.

      INTEGER                                                           &
     & I,J                                                              &
                  ! LOCAL Loop counter (horizontal field index).
     &,K                                                                &
                  ! LOCAL Tile pointer
     &,L                                                                &
                  ! LOCAL Land pointer
     &,N          ! LOCAL Loop counter (tile index).
!== sxy
! local variables in the subroutine sf_evap
      REAL                                                              &
     & DFQW(LAND_PTS)                                                   &
                             ! Increment in GBM moisture flux.
     &,DFTL(LAND_PTS)                                                   &
                             ! Increment in GBM sensible heat flux.
     &,E_TILE_OLD(LAND_PTS,NTILES)                                      &
!                            ! Surface moisture flux before adjustment.
     &,LE_TILE_OLD(LAND_PTS,NTILES)
!                            ! Surf latent heat flux before adjustment.
      REAL                                                              &
     & DIFF_LAT_HTF                                                     &
                             ! Increment in local latent heat flux.
     &,DIFF_SENS_HTF                                                    &
                             ! Increment in local sensible heat flux.
     &,DTSTAR                                                           &
                             ! Increment in local surface temperature.
     &,EDT                   ! Moisture flux x timestep

      INTEGER                                                           &
     & M                    ! Loop counter (soil level index).          
! local variables in the subroutine sf_melt
      REAL                                                              &
     & DFQW_M                                                             &
                            ! Moisture flux increment.
     &,DFTL_M                                                             &
                            ! Sensible heat flux increment.
!     &,DTSTAR                                                           &
                            ! Surface temperature increment.
     &,LCMELT                                                           &
                            ! Temporary in melt calculations.
     &,LSMELT                                                           &
                            ! Temporary in melt calculations.
     &,RHOKH1_PRIME                                                     &
                            ! Modified forward time-weighted
!                           ! transfer coefficient.
     &,SNOW_MAX                                                         &
                            ! Snow available for melting.
     &,TSTARMAX             ! Maximum gridbox mean surface temperature
!                           ! at sea points with ice. 
!jhan{:delete tis diag. sect
! EAK
!      INTEGER, PARAMETER  :: ntest = 0 !  for prints


      REAL ltfs

!     GS_TILE (in CABLE) is returned to GC in the UM (pfv, 25oct13a)
!     REAL                                                              &
!    & GS_TILE(LAND_PTS,NTILES)  ! returned from CABLE 


      ltfs = access_tfs

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_IMPL ',3)
      ENDIF
      ERROR = 0

! DEPENDS ON: im_sf_pt
      CALL IM_SF_PT (                                                   &
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,LAND_PTS                      &
     &,LAND_INDEX,NTILES,TILE_INDEX,TILE_PTS                            &
     &,FLANDG,TILE_FRAC,SNOW_TILE,ICE_FRACT                             &
     &,GAMMA,ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE                        &
     &,RESFT,RHOKPM,RHOKPM_POT,RHOKPM_SICE                              &
     &,RHOKM_U_1,RHOKM_V_1,RHOKH_TILE,RHOKH_SICE                        &
     &,CT_CTQ_1,DQW_1,DTL_1,CQ_CM_U_1,CQ_CM_V_1,DU_1,DV_1               &
     &,FLANDG_U,FLANDG_V                                                &
     &,FQW_1,FTL_1                                                      &
     &,TAUX_1,TAUX_LAND,TAUX_SSI,TAUY_1,TAUY_LAND,TAUY_SSI              &
     &,FQW_TILE,EPOT_TILE,FTL_TILE,FQW_ICE,FTL_ICE,E_SEA,H_SEA          &
     &,L_FLUX_BC,LTIMER                                                 &
! EAK 
     &,l_cable                                                          &
     &)

!-----------------------------------------------------------------------
!! 6.1 Convert FTL to sensible heat flux in Watts per square metre.
!-----------------------------------------------------------------------

!fpp$ Select(CONCUR)
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        FTL_1(I,J) = FTL_1(I,J)*CP
       ENDDO
      ENDDO

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        FTL_ICE(I,J) = CP*FTL_ICE(I,J)
       ENDDO
      ENDDO

      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          FTL_TILE(L,N) = CP*FTL_TILE(L,N)
        ENDDO
      ENDDO
      DO N=1,NTILES
        DO L=1,LAND_PTS
          MELT_TILE(L,N) = 0.
        ENDDO
      ENDDO

!     GS_TILE (in CABLE) is returned to GC in the UM (pfv, 25oct13a)

      if( l_cable ) then
! DEPENDS ON: cable_implicit_driver.o
   call cable_implicit_driver(                      &
                     LS_RAIN, CONV_RAIN, LS_SNOW, CONV_SNOW, dtl_1, dqw_1, &
                     T_SOIL, TSOIL_TILE, SMCL, SMCL_TILE, timestep,                            &
                     SMVCST,STHF, STHF_TILE, STHU,STHU_TILE, snow_tile,                         &
                     SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L, SNOW_MASS3L, &
                     SNOW_RHO3L,SNOW_TMP3L,SNOW_COND, FTL_TILE_CAB,     &
                     FTL_CAB,LE_TILE_CAB,LE_CAB, FTL_1,FTL_TILE,FQW_1,  &
                     FQW_TILE, TSTAR_TILE, TSTAR_TILE_CAB, TSTAR_CAB, &
                     SMCL_CAB,TSOIL_CAB, SURF_HTF_CAB, SURF_HT_FLUX_LAND, &
                     ECAN_TILE,ESOIL_TILE,EI_TILE,RADNET_TILE,TOT_ALB,  &
                     SNAGE_TILE,CANOPY, GS, GC, T1P5M_TILE, Q1P5M_TILE, &
                     CANOPY_GB,FLAND,MELT_TILE, DIM_CS1,DIM_CS2, NPP, NPP_FT, &
                     GPP, GPP_FT, RESP_S, RESP_S_TOT, RESP_S_TILE,     &
                     RESP_P,RESP_P_FT, G_LEAF, TRANSP_TILE, CPOOL_TILE,       &
                     NPOOL_TILE, PPOOL_TILE, GLAI, PHENPHASE, &
                     NPP_FT_ACC, RESP_W_FT_ACC, iday_number &
                      )

   endif

!-----------------------------------------------------------------------
! Diagnose the GBM surface temperature for points with sea-ice
!-----------------------------------------------------------------------

      TSTAR_SIC(:,:,:)= 0.0
      TSTAR_SSI0(:,:) = TSTAR_SSI(:,:)
      SURF_HT_FLUX_SICE(:,:)=0.0

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        IF ( FLANDG(I,J) <  1.0 .and. ICE_FRACT(I,J) >  0. ) THEN
          SURF_HT_FLUX_SICE(I,J) = RADNET_SICE(I,J) - LS*FQW_ICE(I,J) - &
     &                        FTL_ICE(I,J)
          TSTAR_SSI0(I,J) = (1. - ICE_FRACT(I,J))*TSTAR_SEA(I,J) +      &
     &                   ICE_FRACT(I,J)*TI(I,J,1) +                     &
     &                    SURF_HT_FLUX_SICE(I,J) / ASHTF(I,J)
          DO N=1,NICE
            TSTAR_SIC(I,J,N) = ICE_FRACT_NCAT(I,J,N)*TI(I,J,N) +        &
     &                   (ICE_FRACT_NCAT(I,J,N)/ICE_FRACT(I,J))*        &
     &                   (SURF_HT_FLUX_SICE(I,J)/ASHTF(I,J))
          ENDDO
        ENDIF

       ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Optional error check : test for negative top soil layer temperature
!-----------------------------------------------------------------------
      IF (L_NEG_TSTAR) THEN
        DO L=1,LAND_PTS
          IF (T_SOIL(L,1) <  0) THEN
            ERROR = 1
            WRITE(6,*) '*** ERROR DETECTED BY ROUTINE BDY_LAYR ***'
            WRITE(6,*) 'NEGATIVE TEMPERATURE IN TOP SOIL LAYER AT '
            WRITE(6,*) 'LAND POINT ',L
          ENDIF
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
!!   Diagnose the land surface temperature
!-----------------------------------------------------------------------

      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          TSTAR_TILE_OLD(L,N) = TSTAR_TILE(L,N)
        ENDDO
      ENDDO

! EAK
      IF( .NOT. l_cable ) THEN
       DO N=1,NTILES
        DO L=1,LAND_PTS
          TSTAR_TILE(L,N) = T_SOIL(L,1)
          IF (SNOW_TILE(L,N) >  0.)                                     &
     &      TSTAR_TILE(L,N) =  MIN( T_SOIL(L,1), TM )
        ENDDO
       ENDDO

       DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          LAT_HT = LC
          IF (SNOW_TILE(L,N) >  0.) LAT_HT = LS
          TSTAR_TILE(L,N) = T_SOIL(L,1) + ( RADNET_TILE(L,N)            &
     &                          - LAT_HT*FQW_TILE(L,N) - FTL_TILE(L,N)  &
     &                          + (CANHC_TILE(L,N)/TIMESTEP) *          &
     &                            (TSTAR_TILE_OLD(L,N) - T_SOIL(L,1)) ) &
     &                                                 / ASHTF_TILE(L,N)
        ENDDO
       ENDDO
      ENDIF  ! IF( .NOT. l_cable )

!-----------------------------------------------------------------------
!! 7.  Surface evaporation components and updating of surface
!!     temperature (P245, routine SF_EVAP).
!-----------------------------------------------------------------------
! 
      IF(  l_cable ) THEN
       DO N=1,NTILES
        DO L=1,LAND_PTS
!          ECAN_TILE(L,N) = 0.
!          ESOIL_TILE(L,N) = 0.
          ELAKE_TILE(L,N) = 0.
          EI_TILE(L,N) = 0.
        ENDDO
       ENDDO

       DO J=1,ROWS
        DO I=1,ROW_LENGTH
         ECAN(I,J) = 0.
         ESOIL(I,J) = 0.
        ENDDO
       ENDDO

       DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
!          E_TILE_OLD(L,N) = FQW_TILE(L,N)
!          IF (SNOW_TILE(L,N)  >   0.) THEN
!            LE_TILE_OLD(L,N) = (LC + LF)*FQW_TILE(L,N)
!          ELSE
!            LE_TILE_OLD(L,N) = LC*FQW_TILE(L,N)
!          ENDIF
          ! Sublimation from snow-covered land tiles
!          IF (SNOW_TILE(L,N)  >   0.) THEN
!            EI_TILE(L,N) =  FQW_TILE(L,N)
!            EDT = EI_TILE(L,N)*TIMESTEP
!            IF ( EDT  >   SNOW_TILE(L,N) )                              &
!     &        EI_TILE(L,N) = SNOW_TILE(L,N) / TIMESTEP
!            FQW_TILE(L,N) = FQW_TILE(L,N) -  EI_TILE(L,N)
!          ENDIF
          ! Surface evaporation from and condensation onto snow-free land
!          IF ( FQW_TILE(L,N)  >   0.0 ) THEN
!            ECAN_TILE(L,N) = (1. - FLAKE(L,N)) *                        &
!     &                       FRACA(L,N) * FQW_TILE(L,N) / RESFT(L,N)
!            ESOIL_TILE(L,N) = (1. - FLAKE(L,N)) *                       &
!     &                        (1. - FRACA(L,N))*RESFS(L,N)*FQW_TILE(L,N)&
!     &                                                      / RESFT(L,N)
!            ELAKE_TILE(L,N) = FLAKE(L,N)*FQW_TILE(L,N) / RESFT(L,N)
!            EDT = ECAN_TILE(L,N)*TIMESTEP
!            IF ( EDT  >   CANOPY(L,N) ) THEN
!              ESOIL_TILE(L,N) =  (1. - FLAKE(L,N)) *                    &
!     &                           (1. - FRACA(L,N)*CANOPY(L,N)/EDT) *    &
!     &                               RESFS(L,N)*FQW_TILE(L,N)/RESFT(L,N)
!              ECAN_TILE(L,N) = CANOPY(L,N) / TIMESTEP
!            ENDIF
!          ELSEIF (SNOW_TILE(L,N) <= 0.) THEN
!            IF (TSTAR_TILE(L,N) >= TM) THEN
!              ECAN_TILE(L,N) = (1. - FLAKE(L,N))*FQW_TILE(L,N)
!              ELAKE_TILE(L,N) = FLAKE(L,N)*FQW_TILE(L,N)
!            ELSE
!              EI_TILE(L,N) =  FQW_TILE(L,N)
!            ENDIF
!          ENDIF

          ECAN(I,J) = ECAN(I,J) + TILE_FRAC(L,N)*ECAN_TILE(L,N)
          ESOIL(I,J) = ESOIL(I,J) + TILE_FRAC(L,N)*ESOIL_TILE(L,N)          
        ENDDO
       ENDDO   
       EXT = 0. ! MRD
!-----------------------------------------------------------------------    
! Soil evapotranspiration
!-----------------------------------------------------------------------
!       DO L=1,LAND_PTS
!         J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
!         I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
!         EDT = ESOIL(I,J)*TIMESTEP
!         IF ( EDT  >   SMC(L) ) THEN
!           DO N=1,NTILES
!             ESOIL_TILE(L,N) = SMC(L)*ESOIL_TILE(L,N) / EDT
!           ENDDO
!           ESOIL(I,J) = SMC(L) / TIMESTEP
!         ENDIF
!       ENDDO

!       DO M=1,SM_LEVELS
!         DO L=1,LAND_PTS
!           EXT(L,M) = 0.
!         ENDDO
!       ENDDO
!       DO M=1,SM_LEVELS
!         DO N=1,NTILES
!           DO K=1,TILE_PTS(N)
!            L = TILE_INDEX(K,N)
!            EXT(L,M) = EXT(L,M) + TILE_FRAC(L,N)*WT_EXT_TILE(L,M,N)     &
!     &                                          *ESOIL_TILE(L,N)
!           ENDDO
!         ENDDO
!       ENDDO
!-----------------------------------------------------------------------
! Calculate increments to surface heat fluxes, moisture fluxes and
! temperatures
!-----------------------------------------------------------------------
!       DO L=1,LAND_PTS
!         DFTL(L) = 0.
!         DFQW(L) = 0.
!       ENDDO
!
!       DO N=1,NTILES
!         DO K=1,TILE_PTS(N)
!           L = TILE_INDEX(K,N)
!           DIFF_LAT_HTF = (LC + LF)*EI_TILE(L,N) + LC*ECAN_TILE(L,N)     &
!     &                     + LC*ESOIL_TILE(L,N) + LC*ELAKE_TILE(L,N)     &
!     &                    - LE_TILE_OLD(L,N)
!           DIFF_SENS_HTF = - DIFF_LAT_HTF /                              &
!     &                        ( 1. + ASHTF_TILE(L,N)/(CP*RHOKH_TILE(L,N)) )
!           FTL_TILE(L,N) = FTL_TILE(L,N) + DIFF_SENS_HTF
!           DTSTAR = - (DIFF_LAT_HTF + DIFF_SENS_HTF) / ASHTF_TILE(L,N)
!           TSTAR_TILE(L,N) = TSTAR_TILE(L,N) + DTSTAR
!           DFTL(L) = DFTL(L) + TILE_FRAC(L,N)*DIFF_SENS_HTF
!           DFQW(L) = DFQW(L) + TILE_FRAC(L,N)*( ECAN_TILE(L,N) +         &
!     &                   ESOIL_TILE(L,N) + EI_TILE(L,N) + ELAKE_TILE(L,N)&
!     &                   - E_TILE_OLD(L,N) )
!         ENDDO
!       ENDDO
!!-----------------------------------------------------------------------
! Update level 1 temperature and humidity and GBM heat and moisture
! fluxes due to limited moisture availability
!-----------------------------------------------------------------------
!       DO L=1,LAND_PTS
!         J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
!         I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
!         FTL_1(I,J) = FTL_1(I,J) + FLAND(L)*DFTL(L)
!         FQW_1(I,J) = FQW_1(I,J) + FLAND(L)*DFQW(L)
!       ENDDO
 
      ELSE ! ( l_cable )

! DEPENDS ON: sf_evap
      CALL SF_EVAP (                                                    &
     & ROW_LENGTH,ROWS,LAND_PTS,NTILES,                                 &
     & LAND_INDEX,TILE_INDEX,TILE_PTS,SM_LEVELS,LTIMER,FLAND,           &
     & ASHTF_TILE,CANOPY,DTRDZ_CHARNEY_GRID_1,FLAKE,FRACA,SNOW_TILE,    &
     & RESFS,RESFT,RHOKH_TILE,TILE_FRAC,SMC,WT_EXT_TILE,TIMESTEP,       &
     & FQW_1,FQW_TILE,FTL_1,FTL_TILE,TSTAR_TILE,                        &
     & ECAN,ECAN_TILE,ELAKE_TILE,ESOIL,ESOIL_TILE,EI_TILE,EXT           &
     & )
! EAK
      ENDIF     ! ( l_cable ) THEN

!-----------------------------------------------------------------------
!!     Surface melting of sea-ice and snow on land tiles.
!-----------------------------------------------------------------------
! 
     IF(  l_cable ) THEN
       DO J=1,ROWS
        DO I=1,ROW_LENGTH
         SICE_MELT0(I,J) = 0.0
         DO N=1,NICE
           SICE_MELT(I,J,N) = 0.0
           IF (SIMLT) SICE_MLT_HTF(I,J,N) = 0.0
         ENDDO
         SNOMLT_SURF_HTF(I,J) = 0.0
         SNOWMELT(I,J) = 0.0
         EI_LAND(I,J) = 0.0
         EI_SICE(I,J) = 0.0
        ENDDO
       ENDDO

!       DO N=1,NTILES
!         DO L=1,LAND_PTS
!           MELT_TILE(L,N) = 0.
!         ENDDO
!       ENDDO
!       SNOMLT_SURF_HTF = 0. ! MRD
!-----------------------------------------------------------------------
!  Melt snow on land tiles if TSTAR_TILE is greater than TM.
!-----------------------------------------------------------------------
!       DO N=1,NTILES
!         DO K=1,TILE_PTS(N)
!           L = TILE_INDEX(K,N)
!           J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
!           I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
!           SNOW_MAX = MAX( 0.0, SNOW_TILE(L,N) - EI_TILE(L,N)*TIMESTEP )
!           IF ( SNOW_MAX >  0.0 .AND. TSTAR_TILE(L,N) >  TM ) THEN
!             LCMELT = (CP + LC*ALPHA1(L,N)*RESFT(L,N))*RHOKH_TILE(L,N)      &
!     &               + ASHTF_TILE(L,N)
!             LSMELT = LCMELT + LF*ALPHA1(L,N)*RHOKH_TILE(L,N)
!             DTSTAR = - MIN( TSTAR_TILE(L,N) - TM ,                      &
!     &                      LF*SNOW_MAX / (LCMELT*TIMESTEP) )
!             MELT_TILE(L,N) = - LSMELT*DTSTAR / LF
!             DFTL_M = CP*RHOKH_TILE(L,N)*DTSTAR
!             DFQW_M = ALPHA1(L,N)*RESFT(L,N)*RHOKH_TILE(L,N)*DTSTAR
!             FTL_TILE(L,N) = FTL_TILE(L,N) + DFTL_M
!             EI_TILE(L,N) = EI_TILE(L,N) + DFQW_M
!             TSTAR_TILE(L,N) = TSTAR_TILE(L,N) + DTSTAR       
!!-----------------------------------------------------------------------
!!  Update gridbox-mean quantities
!!-----------------------------------------------------------------------
!             DFTL_M = TILE_FRAC(L,N)*DFTL_M
!             DFQW_M = TILE_FRAC(L,N)*DFQW_M
!             FTL_1(I,J) = FTL_1(I,J) + FLANDG(I,J)*DFTL_M
!             FQW_1(I,J) = FQW_1(I,J) + FLANDG(I,J)*DFQW_M
!           ENDIF
!           EI_LAND(I,J) = EI_LAND(I,J) + TILE_FRAC(L,N)*EI_TILE(L,N)
!         ENDDO
!       ENDDO
!-----------------------------------------------------------------------
!  Increment snow by sublimation and melt
!-----------------------------------------------------------------------
       DO N=1,NTILES
         DO K=1,TILE_PTS(N)
           L = TILE_INDEX(K,N)
           J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
           I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
!           SNOW_TILE(L,N) = SNOW_TILE(L,N) -                             &
!     &                     (EI_TILE(L,N) + MELT_TILE(L,N))*TIMESTEP
           SNOWMELT(I,J) = SNOWMELT(I,J) + TILE_FRAC(L,N)*MELT_TILE(L,N)
         ENDDO
       ENDDO
!       IF (SMLT) THEN
       DO J=1,ROWS
         DO I=1,ROW_LENGTH
          SNOMLT_SURF_HTF(I,J) = LF*SNOWMELT(I,J)
         ENDDO
       ENDDO
!       ENDIF

       IF(NICE  ==  1)THEN
       DO J=1,ROWS
        DO I=1,ROW_LENGTH
          DTSTAR=0.0
         IF ( FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0.0 ) THEN
!-----------------------------------------------------------------------
!   Melt sea-ice if TSTAR > TSTARMAX
!-----------------------------------------------------------------------
           EI_SICE(I,J) = FQW_ICE(I,J)
           TSTARMAX = ICE_FRACT(I,J)*TM                                  &
     &         + (1.0 - ICE_FRACT(I,J))*TSTAR_SEA(I,J)
           IF ( TSTAR_SSI0(I,J)  >   TSTARMAX ) THEN
             RHOKH1_PRIME = 1. / ( 1. / RHOKH_SICE(I,J)                &
     &                       + ICE_FRACT(I,J)*GAMMA*DTRDZ_CHARNEY_GRID_1(I,J) )
             DTSTAR = TSTARMAX - TSTAR_SSI0(I,J)
             LSMELT = (CP + (LC + LF)*ALPHA1_SICE(I,J))*RHOKH1_PRIME     &
     &                                                 + ASHTF(I,J)
             DFTL_M = CP * RHOKH1_PRIME * DTSTAR
             DFQW_M = ALPHA1_SICE(I,J) * RHOKH1_PRIME * DTSTAR
             TSTAR_SSI0(I,J) = TSTARMAX
             SICE_MELT0(I,J) = - LSMELT * DTSTAR
             IF (SIMLT) SICE_MLT_HTF(I,J,1) = SICE_MELT0(I,J)
             FTL_1(I,J) = FTL_1(I,J) + (1.0-FLANDG(I,J))*DFTL_M
             FQW_1(I,J) = FQW_1(I,J) + (1.0-FLANDG(I,J))*DFQW_M
             EI_SICE(I,J) = EI_SICE(I,J) + DFQW_M
             FTL_ICE(I,J) = FTL_ICE(I,J) + DFTL_M
             FQW_ICE(I,J) = FQW_ICE(I,J) + DFQW_M

           ENDIF

         ENDIF
        ENDDO
       ENDDO

       ELSE

        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          IF ( FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0.0 ) THEN   
!-----------------------------------------------------------------------
!   Melt sea-ice if TSTAR > TSTARMAX
!-----------------------------------------------------------------------
            EI_SICE(I,J) = FQW_ICE(I,J)
            DTSTAR=0.0
            DO N=1,NICE
              IF (ICE_FRACT_NCAT(I,J,N) >  0.0) THEN
                TSTARMAX = ICE_FRACT_NCAT(I,J,N)*TM
                IF (TSTAR_SIC(I,J,N) >  TSTARMAX) THEN
                  RHOKH1_PRIME = 1. / ( 1. / RHOKH_SICE(I,J)           &
     &                 + ICE_FRACT_NCAT(I,J,N)*GAMMA*DTRDZ_CHARNEY_GRID_1(I,J) )
                  DTSTAR = TSTARMAX - TSTAR_SIC(I,J,N)
                  LSMELT = (CP + (LC + LF)*ALPHA1_SICE(I,J))*RHOKH1_PRIME&
     &                                                 + ASHTF(I,J)
                  DFTL_M = CP * RHOKH1_PRIME * DTSTAR
                  DFQW_M = ALPHA1_SICE(I,J) * RHOKH1_PRIME * DTSTAR
                  TSTAR_SIC(I,J,N) = TSTARMAX
                  SICE_MELT(I,J,N) = - LSMELT * DTSTAR
                  IF (SIMLT) SICE_MLT_HTF(I,J,N) = SICE_MELT(I,J,N)
                  FTL_1(I,J)   = FTL_1(I,J) + (1.0-FLANDG(I,J))*DFTL_M
                  FQW_1(I,J)   = FQW_1(I,J) + (1.0-FLANDG(I,J))*DFQW_M
                  EI_SICE(I,J) = EI_SICE(I,J) + DFQW_M
                  FTL_ICE(I,J) = FTL_ICE(I,J) + DFTL_M
                  FQW_ICE(I,J) = FQW_ICE(I,J) + DFQW_M

                ENDIF
              ENDIF
            ENDDO
          ENDIF
         ENDDO
        ENDDO
       ENDIF

      ELSE !( l_cable )

! DEPENDS ON: sf_melt
      CALL SF_MELT (                                                    &
     & ROW_LENGTH,ROWS,LAND_PTS,NTILES,NICE,LAND_INDEX,                 &
     & TILE_INDEX,TILE_PTS,LTIMER,SIMLT,SMLT,FLANDG,                    &
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,DTRDZ_CHARNEY_GRID_1,        &
     & ICE_FRACT,RESFT,RHOKH_TILE,RHOKH_SICE,TILE_FRAC,TIMESTEP,        &
     & GAMMA,                                                           &
     & EI_TILE,FQW_1,FQW_ICE,FTL_1,FTL_ICE,FTL_TILE,                    &
     & TSTAR_SEA,TSTAR_SIC,TSTAR_TILE,SNOW_TILE,                        &
     & EI_LAND,EI_SICE,SICE_MELT,ICE_FRACT_NCAT,                        &
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,MELT_TILE,                 &
     & TSTAR_SSI0,SICE_MELT0)

      ENDIF ! ( l_cable ) THEN


      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        QIM_1(I,J)=QW_1(I,J) + DQW_1(I,J)-CT_CTQ_1(I,J)*FQW_1(I,J)
        TIM_1(I,J)=TL_1(I,J) + DTL_1(I,J)-CT_CTQ_1(I,J)*FTL_1(I,J)/CP
        IF ( FLANDG(I,J) <  1.0 .and. ICE_FRACT(I,J) >  0. ) THEN
          IF(NICE  ==  1)THEN
            TSTAR_SSI(I,J)=TSTAR_SSI0(I,J)
          ELSE
            TSTAR_SSI(I,J)=(1.-ICE_FRACT(I,J))*TSTAR_SEA(I,J)
            DO N=1,NICE
              TSTAR_SSI(I,J)=TSTAR_SSI(I,J)+TSTAR_SIC(I,J,N)
            ENDDO
          ENDIF
        ENDIF
       ENDDO
      ENDDO

!jhan{:delete tis diag. sect
!      IF(ntest>0) THEN
!      print *,'l=100,i=4,j=2 mp '
!      l=100
!      i=4
!      j=2
!      N=9
!      print 21,l,SURF_RADFLUX(I,J),RADNET_TILE(l,N),DQW_1(I,J), &
!                                         DTL_1(I,J),CT_CTQ_1(I,J)
!      print 22,l,SNOW_TILE(l,N),TL_1(i,j),T_SOIL(l,1:4)
!      print 23,l,FQW_1(i,j),FQW_TILE(l,N)
!      print 24,l,FTL_1(i,j),FTL_TILE(l,N),FLAND(l),TILE_FRAC(l,N)
!      print 25,l,Z0H_TILE(l,N),Z0M_TILE(l,N)
!      print 26,l,TSTAR(i,j),TSTAR_TILE(l,N)
!21    format(x,'sIRADNET',i4,'z',2f7.1,x,2f8.4,2x,f10.5)
!22    format(x,'sIRHOKH',i4,'z',f8.1,x,f5.1,x,4f6.1)
!23    format(x,'sIFQW_1',i4,'z',f8.5,x,f8.5,x,f5.3,x,f4.0,x,f4.0)
!24    format(x,'sIFTL_1',i4,'z',f8.2,x,f8.2,x,f5.3,x,f6.3,x,f6.2)
!25    format(x,'sIZ0M',i4,'z',f8.5,x,f8.5,x,f5.2)
!26    format(x,'sITSTAR',i4,'z',f6.1,x,f6.1,x,f7.5,x,f5.0)
!
!      print *,'l=1000,i=13,j=29'
!      l=1000
!      i=13
!      j=29
!      l=2132
!      i=28
!      j=57
!      print 31,l,SURF_RADFLUX(I,J),RADNET_TILE(l,1:5),RADNET_TILE(l,8), &
!                  DQW_1(I,J), DTL_1(I,J),CT_CTQ_1(I,J)
!      print 32,l,SNOW_TILE(l,1:5),SNOW_TILE(l,8),TL_1(i,j),T_SOIL(l,1:4)
!      print 33,l,FQW_1(i,j),FQW_TILE(l,1:5),FQW_TILE(l,8)
!      print 34,l,FTL_1(i,j),FTL_TILE(l,1:5),FTL_TILE(l,8), &
!               FLAND(l),TILE_FRAC(l,1:5),TILE_FRAC(l,8)
!      print 35,l,Z0H_TILE(l,1:5),Z0H_TILE(l,8), &
!               Z0M_TILE(l,1:5),Z0M_TILE(l,8)
!      print 36,l,TSTAR_LAND(i,j),TSTAR_TILE(l,1:5),TSTAR_TILE(l,8)
!31    format(x,'sIRADNET',i4,'z',7f7.1,x,2f8.4,2x,f10.5)
!32    format(x,'sIRHOKH',i4,'z',6f7.1,x,f5.1,x,4f6.1)
!33    format(x,'sIFQW_1',i4,'z',f9.6,x,6f6.4,x,6f5.3,x,6f4.0,x,6f4.0)
!34    format(x,'sIFTL_1',i4,'z',f8.3,x,6f6.1,x,f5.3,x,6f6.3,x,6f4.2)
!35    format(x,'sIZ0M',i4,'z',6f7.5,x,6f7.5,x,6f4.2)
!36    format(x,'sITSTAR',i4,'z',f6.1,x,6f6.1,x,6f7.5,x,6f6.0)
!!     l=1500 has i=94,j=44
!      print *,'l=1500,i=94,j=44'
!      l=1500
!      i=94
!      j=44
!      l=2080
!      i=35
!      j=56 
!
!      print 31,l,SURF_RADFLUX(I,J),RADNET_TILE(l,1:5),RADNET_TILE(l,8),DQW_1(I,J), &
!               DTL_1(I,J),CT_CTQ_1(I,J)
!      print 32,l,SNOW_TILE(l,1:5),SNOW_TILE(l,8),TL_1(i,j),T_SOIL(l,1:4)
!      print 33,l,FQW_1(i,j),FQW_TILE(l,1:5),FQW_TILE(l,8)
!      print 34,l,FTL_1(i,j),FTL_TILE(l,1:5),FTL_TILE(l,8), &
!               FLAND(l),TILE_FRAC(l,1:5),TILE_FRAC(l,8)
!      print 35,l,Z0H_TILE(l,1:5),Z0H_TILE(l,8), &
!               Z0M_TILE(l,1:5),Z0M_TILE(l,8)
!      print 36,l,TSTAR_LAND(i,j),TSTAR_TILE(l,1:5),TSTAR_TILE(l,8)
!      ENDIF
!jhan}
!-----------------------------------------------------------------------
!!     Specific humidity and temperature at 1.5 metres.
!-----------------------------------------------------------------------
! DEPENDS ON: screen_tq
      CALL SCREEN_TQ (                                                  &
     & ROW_LENGTH,ROWS,LAND_PTS,NTILES,                                 &
     & LAND_INDEX,TILE_INDEX,TILE_PTS,FLANDG,                           &
     & SQ1P5,ST1P5,CHR1P5M,CHR1P5M_SICE,PSTAR,QIM_1,RESFT,              &
     & TILE_FRAC,TIM_1,TSTAR_SSI,TSTAR_TILE,                            &
     & Z0HSSI,Z0H_TILE,Z0MSSI,Z0M_TILE,Z1,                              &
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE,                               &
     & lq_mix_bl,l_cable                                                &
     & )

!-----------------------------------------------------------------------
!!     Gridbox-mean surface temperature and net surface heat fluxes
!-----------------------------------------------------------------------
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        IF( .NOT. l_cable ) SURF_HT_FLUX_LAND(I,J) = 0.
        SURF_HT_FLUX_SICE(I,J) = 0.
        DO N=1,NICE
          SURF_HT_FLUX_SICE_NCAT(I,J,N)=0.
        ENDDO
        TSTAR_RAD4(I,J) = 0.
        IF (FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0.) THEN
          TSTAR_SICE(I,J) = (TSTAR_SSI(I,J) -                           &
     &              (1.-ICE_FRACT(I,J))*TSTAR_SEA(I,J))/ICE_FRACT(I,J)
          TSTAR_RAD4(I,J) = (1.0-FLANDG(I,J))                           &
     &              *ICE_FRACT(I,J)*TSTAR_SICE(I,J)**4
          RADNET_SICE(I,J) = SURF_RADFLUX(I,J) -                            &
     &                     ICE_FRACT(I,J)*SBCON*TSTAR_SICE(I,J)**4
          IF(NICE  ==  1)THEN
            SURF_HT_FLUX_SICE(I,J) = RADNET_SICE(I,J) - LS*FQW_ICE(I,J)-&
     &                        FTL_ICE(I,J) - SICE_MELT0(I,J)
            SURF_HT_FLUX_SICE_NCAT(I,J,1) = SURF_HT_FLUX_SICE(I,J)
          ELSE
            DO N=1,NICE
              SURF_HT_FLUX_SICE_NCAT(I,J,N) = (ICE_FRACT_NCAT(I,J,N) /  &
     &                        ICE_FRACT(I,J)) * (RADNET_SICE(I,J) -     &
     &                        LS*FQW_ICE(I,J) -                         &
     &                        FTL_ICE(I,J)) - SICE_MELT(I,J,N)          &
     &                        + 4.0*SBCON*(TSTAR_SICE(I,J)**3.0) *      &
     &                        (ICE_FRACT_NCAT(I,J,N)*TSTAR_SICE(I,J)    &
     &                        - TSTAR_SIC(I,J,N))
              SURF_HT_FLUX_SICE(I,J) = SURF_HT_FLUX_SICE(I,J) +         &
     &                               SURF_HT_FLUX_SICE_NCAT(I,J,N)
            ENDDO
          ENDIF

        ENDIF
       ENDDO
      ENDDO

      DO L=1,LAND_PTS
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        TSTAR_LAND(I,J) = 0.
      ENDDO

!-----------------------------------------------------------------------
!! The following calculation on land tiles is repeated for 
!! LRAD_EMIS_LAND_GEN equals false for bit reproducibility.
!! It would be tidier to move the if statement within the do loop
!! but this affects vectorisation.
!-----------------------------------------------------------------------
!     calculate net radiation for land points (in surf_radflux); EAK
      DO N=1,NTILES
      DO K=1,TILE_PTS(N)
        L = TILE_INDEX(K,N)
        J =(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        SURF_RADFLUX(I,J) = 0.0 
      ENDDO
      ENDDO


      If (LRAD_EMIS_LAND_GEN) Then

! ANTHROP_HEAT(N) takes value 0.0 except when N=6 (urban) and
! L_ANTHROP_HEAT_SRC is true

        DO N=1,NTILES
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
            IF( .NOT. l_cable ) RADNET_TILE(L,N) = SW_TILE(L,N) +       &
     &                        RAD_EMIS_LAND_GEN*        &
     &                        (LW_DOWN(I,J) - SBCON*TSTAR_TILE(L,N)**4)
            LE_TILE(L,N) = LC*ECAN_TILE(L,N) + LC*ESOIL_TILE(L,N) +     &
     &                     LC*ELAKE_TILE(L,N) + LS*EI_TILE(L,N)
          IF( .NOT. l_cable ) SURF_HT_FLUX_LAND(I,J) =                  &
     &                        SURF_HT_FLUX_LAND(I,J)                    &
     &                      + TILE_FRAC(L,N) *                          &
     &                      ( RADNET_TILE(L,N) + ANTHROP_HEAT(N) -      &
     &                        FTL_TILE(L,N) -                           &
     &                        LE_TILE(L,N) - LF*MELT_TILE(L,N) -        &
     &                       (CANHC_TILE(L,N)/TIMESTEP) *               &
     &                       (TSTAR_TILE(L,N) - TSTAR_TILE_OLD(L,N)) )
            TSTAR_LAND(I,J) = TSTAR_LAND(I,J)                           &
     &                 + TILE_FRAC(L,N)*TSTAR_TILE(L,N)
            TSTAR_RAD4(I,J) = TSTAR_RAD4(I,J) + FLANDG(I,J)*            &
     &                      TILE_FRAC(L,N)*TSTAR_TILE(L,N)**4
            SURF_RADFLUX(I,J) =  SURF_RADFLUX(I,J) + FLANDG(I,J)*       &
                                 TILE_FRAC(L,N)*RADNET_TILE(L,N)
          ENDDO
        ENDDO

      Else

        DO N=1,NTILES
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
            IF( .NOT. l_cable ) RADNET_TILE(L,N) = SW_TILE(L,N) +       &
     &                        LW_DOWN(I,J) - SBCON*TSTAR_TILE(L,N)**4
            LE_TILE(L,N) = LC*ECAN_TILE(L,N) + LC*ESOIL_TILE(L,N) +     &
     &                     LC*ELAKE_TILE(L,N) + LS*EI_TILE(L,N)
          IF( .NOT. l_cable ) SURF_HT_FLUX_LAND(I,J) =                  &
                              SURF_HT_FLUX_LAND(I,J)                    &
     &                      + TILE_FRAC(L,N) *                          &
     &                      ( RADNET_TILE(L,N) + ANTHROP_HEAT(N) -      &
     &                        FTL_TILE(L,N) -                           &
     &                        LE_TILE(L,N) - LF*MELT_TILE(L,N) -        &
     &                       (CANHC_TILE(L,N)/TIMESTEP) *               &
     &                       (TSTAR_TILE(L,N) - TSTAR_TILE_OLD(L,N)) )
            TSTAR_LAND(I,J) = TSTAR_LAND(I,J)                           &
     &                 + TILE_FRAC(L,N)*TSTAR_TILE(L,N)
            TSTAR_RAD4(I,J) = TSTAR_RAD4(I,J) + FLANDG(I,J)*            &
     &                      TILE_FRAC(L,N)*TSTAR_TILE(L,N)**4
            SURF_RADFLUX(I,J) =  SURF_RADFLUX(I,J) + FLANDG(I,J)*       &
                                 TILE_FRAC(L,N)*RADNET_TILE(L,N)

          ENDDO
        ENDDO
      Endif

! TOA outward LW radiation after boundary layer
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        OLR(I,J) = OLR(I,J) + SBCON*TSTAR_RAD4(I,J)
       ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Optional error check : test for negative surface temperature
!-----------------------------------------------------------------------
      IF (L_NEG_TSTAR) THEN
        DO L=1,LAND_PTS
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          IF (TSTAR_LAND(I,J) <  0) THEN
            ERROR = 1
            WRITE(6,*) '*** ERROR DETECTED BY ROUTINE BDY_LAYR ***'
            WRITE(6,*) 'NEGATIVE SURFACE TEMPERATURE AT LAND POINT ',L
          ENDIF
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Update sea-ice surface layer temperature.
!-----------------------------------------------------------------------

! DEPENDS ON: sice_htf
      CALL SICE_HTF(                                                    &
     & ROW_LENGTH,ROWS,FLANDG,SIMLT,NICE,                               &
     & DI_NCAT,ICE_FRACT,ICE_FRACT_NCAT,SURF_HT_FLUX_SICE_NCAT,         &
     & TSTAR_SICE,TIMESTEP,                                             &
     & TI,SICE_MLT_HTF,SEA_ICE_HTF,L_SICE_HEATFLUX,                     &
     & LTIMER)

!----------------------------------------------------------------------
!! 8.1 Update U_V.
!----------------------------------------------------------------------

! U component of 10m wind
      IF (SU10)THEN
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          U10M(I,J) = (U_1(I,J) + (DU_1(I,J) -                          &
     &                   CQ_CM_U_1(I,J)*TAUX_1(I,J)) -                  &
     &                      U_0(I,J))*CDR10M_U(I,J) + U_0(I,J)
         ENDDO
        ENDDO
      ENDIF

! V component of 10m wind
      IF (SV10)THEN
        DO J=1,N_ROWS
         DO I=1,ROW_LENGTH
          V10M(I,J) = (V_1(I,J) + (DV_1(I,J) -                          &
     &                   CQ_CM_V_1(I,J)*TAUY_1(I,J)) -                  &
     &                      V_0(I,J))*CDR10M_V(I,J) + V_0(I,J)
         ENDDO
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
!! 9.  Calculate surface latent heat flux.
!-----------------------------------------------------------------------

      IF (SLH) THEN
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          LATENT_HEAT(I,J) = LC*FQW_1(I,J)                              &
     &     + LF*(FLANDG(I,J)*EI_LAND(I,J)+(1.-FLANDG(I,J))*EI_SICE(I,J))
         ENDDO
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! 10 Set RHOKH, the coefficients required for tracer mixing.
!    Required 5B and after due to change in contents of RHOKH in rest
!    of routine.
!-----------------------------------------------------------------------

      DO J=1,ROWS
        DO I=1,ROW_LENGTH
          RHOKH_MIX(I,J) = RHOKH(I,J)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! When NICE=1 the pointers for Ti_GB and TI are identical and TI=TI_GB
! Thus we need to ensure that neither is reassigned until the end of the
! summing of the gridbox mean. This does not occur when NICE  /=  1 as
! then the pointers are not equivalent.
!-----------------------------------------------------------------------

! At this stage TSTAR from cable is composed from new dt vegetation temp.
! and new soil temp.
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            TSTAR(I,J)=FLANDG(I,J)*TSTAR_LAND(I,J)                      &
     &        +(1.-FLANDG(I,J))*TSTAR_SSI(I,J)
            EI(I,J)=FLANDG(I,J)*EI_LAND(I,J)                            &
     &        +(1.-FLANDG(I,J))*EI_SICE(I,J)

            SURF_HT_FLUX(I,J)=FLANDG(I,J)*SURF_HT_FLUX_LAND(I,J)        &
     &        +(1.-FLANDG(I,J))*SURF_HT_FLUX_SICE(I,J)

          IF (NICE  >   1)THEN
             IF (FLANDG(I,J) == 1.0) THEN
                TI_GB(I,J) = RMDI
             ELSE IF (ICE_FRACT(I,J) <= 0.0) THEN
                if (ocn_sss) then
                    ltfs = ZeroDegC - 0.054 * auscom_salinity(I,J)
                end if
                TI_GB(I,J) = LTFS
             ELSE
                TI_GB(I,J) = 0.0
                DO N=1,NICE
                  TI_GB(I,J)=TI_GB(I,J)+ICE_FRACT_NCAT(I,J,N)           &
     &                  /ICE_FRACT(I,J)*TI(I,J,N)
                ENDDO
             ENDIF
          ENDIF
        ENDDO
      ENDDO

!jhan{:delete tis diag. sect
!      IF(ntest>0) THEN
!      print 29,ftl_1(4,2), ftl_tile(100,9), &
!                fqw_1(4,2),fqw_tile(100,9),CT_CTQ_1(4,2), &
!                TSTAR(4,2),TSTAR_LAND(4,2),TSTAR_TILE(100,9)
!29    format('sfimptheend1',2f8.3,2x,2f10.6,x,f10.4,x,3f6.1)     
!      print 37,ftl_1(28,57), ftl_tile(2132,1:5), &
!               fqw_1(28,57), fqw_tile(2132,1:5),CT_CTQ_1(28,57), &
!      TSTAR(28,57),TSTAR_LAND(28,57),TSTAR_TILE(2132,1:5),TSTAR_TILE(2132,8)
!      print 37,ftl_1(35,56), ftl_tile(2080,1:5), &
!               fqw_1(35,56), fqw_tile(2080,1:5),CT_CTQ_1(35,56), &
!      TSTAR(35,56),TSTAR_LAND(35,56),TSTAR_TILE(2080,1:5),TSTAR_TILE(2080,8)
!37    format('sfimptheend2',6f6.1,x,6f7.6,x,f6.1,8f6.1)     
!!      print 37,ftl_1(13,29), ftl_tile(1000,1:5), &
!!               fqw_1(13,29), fqw_tile(1000,1:5),CT_CTQ_1(13,29), &
!!               TSTAR(13,29),TSTAR_LAND(13,29),TSTAR_TILE(1000,1:5)
!!37    format(1x,'sf_imp theend',6f8.3,2x,6f10.6,x,f10.4,x,7f6.1)     
!!      print 37,ftl_1(94,44), ftl_tile(1500,1:5), &
!!               fqw_1(94,44), fqw_tile(1500,1:5),CT_CTQ_1(94,44), &
!!               TSTAR(94,44),TSTAR_LAND(94,44),TSTAR_TILE(1500,1:5)
!!      print 38,ftl_1(28,57), ftl_tile(2132,1:5), &
!!               fqw_1(28,57), fqw_tile(2132,1:5),CT_CTQ_1(28,57), &
!!               TSTAR(28,57),TSTAR_LAND(28,57),TSTAR_TILE(2132,1:5)
!!38    format(1x,'sf_imp theend 2132',6f8.3,2x,6f10.6,x,f10.4,x,7f6.1)     
!      ENDIF
!!      IF(ntest>0) THEN
!      print *,'l=100,i=4,j=2'
!      l=100
!      i=4
!      j=2
!      N=17
!      LAT_HT = LS
!      print 51,l,SURF_RADFLUX(I,J),RADNET_TILE(l,N),DQW_1(I,J), &
!                                         DTL_1(I,J),CT_CTQ_1(I,J)
!      print 52,l,SNOW_TILE(l,N),QW_1(i,j),TL_1(i,j),T_SOIL(l,1:6)
!      print 53,l,LAT_HT*FQW_1(i,j),LAT_HT*FQW_TILE(l,N)
!      print 54,l,FTL_1(i,j),FTL_TILE(l,N),FLAND(l),TILE_FRAC(l,N)
!      print 55,l,Z0H_TILE(l,N),Z0M_TILE(l,N)
!      print 56,l,TSTAR(i,j),TSTAR_TILE(l,N),SURF_HT_FLUX_LAND(I,J)
!51    format(x,'sIRADNET',i4,'z',2f7.1,x,2f8.4,2x,f10.4)
!52    format(x,'sIRHOKH',i4,'z',f8.1,x,f10.7,x,f5.1,x,6f6.1)
!53    format(x,'sIFQW_1',i4,'z',f8.5,x,f8.5,x,f5.3,x,f4.0,x,f4.0)
!54    format(x,'sIFTL_1',i4,'z',f8.2,x,f8.2,x,f5.3,x,f6.3,x,f6.2)
!55    format(x,'sIZ0M',i4,'z',f10.7,x,f10.7,x,f5.2)
!56    format(x,'sITSTAR',i4,'z',f6.1,x,f6.1,x,f7.1)
!
!      print *,'l=1000,i=13,j=29'
!      l=1000
!      i=13
!      j=29
!      l=2132
!      i=28
!      j=57
!      ! only tiles 3 and 6 are non  zero for l=1000
!      LAT_HT1 = LC
!      IF (SNOW_TILE(L,3) >  0.) LAT_HT1 = LS
!      LAT_HT2 = LC
!      IF (SNOW_TILE(L,6) >  0.) LAT_HT2 = LS
!      print 61,l,SURF_RADFLUX(I,J),RADNET_TILE(l,1:5),RADNET_TILE(l,6), &
!                  DQW_1(I,J), DTL_1(I,J),CT_CTQ_1(I,J)
!      print 62,l,SNOW_TILE(l,1:5),SNOW_TILE(l,6),QW_1(i,j),TL_1(i,j),T_SOIL(l,1:6)
!      print 63,l,LAT_HT1*FQW_1(i,j),LAT_HT1*FQW_TILE(l,1:5),LAT_HT2*FQW_TILE(l,6)
!      print 64,l,FTL_1(i,j),FTL_TILE(l,1:5),FTL_TILE(l,6), &
!               FLAND(l),TILE_FRAC(l,1:6)
!!               FLAND(l),TILE_FRAC(l,1:17)
!      print 65,l,Z0H_TILE(l,1:5),Z0H_TILE(l,6), &
!               Z0M_TILE(l,1:5),Z0M_TILE(l,6)
!      print 66,l,TSTAR_LAND(i,j),TSTAR_TILE(l,1:5),TSTAR_TILE(l,6), &
!                SURF_HT_FLUX_LAND(I,J)
!61    format(x,'sIRADNET',i4,'z',7f7.1,x,2f8.4,2x,f10.5)
!62    format(x,'sIRHOKH',i4,'z',6f7.1,x,f10.7,x,f5.1,x,6f6.1)
!63    format(x,'sIFQW_1',i4,'z',f9.6,x,6f6.4,x,6f5.3,x,6f4.0,x,6f4.0)
!64    format(x,'sIFTL_1',i4,'z',f8.3,x,6f6.1,x,f5.3,x,6f6.3,x,17f4.2)
!65    format(x,'sIZ0M',i4,'z',6f7.5,x,6f7.5,x,6f4.2)
!66    format(x,'sITSTAR',i4,'z',f6.1,x,6f6.1,x,6f6.0)
!!     l=1500 has i=94,j=44
!      print *,'l=1500,i=94,j=44'
!      l=1500
!      i=94
!      j=44
!      l=2080
!      i=35
!      j=56 
!!     here tiles 1, 4 and 9 are non zero 
!      LAT_HT1 = LC
!      IF (SNOW_TILE(L,1) >  0.) LAT_HT1 = LS
!      LAT_HT2 = LC
!      IF (SNOW_TILE(L,4) >  0.) LAT_HT2 = LS
!      LAT_HT3 = LC
!      IF (SNOW_TILE(L,9) >  0.) LAT_HT3 = LS
!
!      print 61,l,SURF_RADFLUX(I,J),RADNET_TILE(l,1:5),RADNET_TILE(l,9),DQW_1(I,J), &
!               DTL_1(I,J),CT_CTQ_1(I,J)
!      print 62,l,SNOW_TILE(l,1:5),SNOW_TILE(l,9),QW_1(i,j),TL_1(i,j),T_SOIL(l,1:6)
!      print 63,l,LAT_HT1*FQW_1(i,j),LAT_HT2*FQW_TILE(l,1:5),LAT_HT3*FQW_TILE(l,9)
!      print 64,l,FTL_1(i,j),FTL_TILE(l,1:5),FTL_TILE(l,9), &
!               FLAND(l),TILE_FRAC(l,1:9)
!      print 65,l,Z0H_TILE(l,1:5),Z0H_TILE(l,9), &
!               Z0M_TILE(l,1:5),Z0M_TILE(l,9)
!      print 66,l,TSTAR_LAND(i,j),TSTAR_TILE(l,1:5),TSTAR_TILE(l,9), &
!               SURF_HT_FLUX_LAND(I,J)
!      ENDIF
!jhan}

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_IMPL ',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_IMPL
