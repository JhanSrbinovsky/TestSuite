#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE SF_IMPL2-----------------------------------------------
!!!
!!!  Purpose: Calculate implicit correction to surface fluxes of heat,
!!!           moisture and momentum to be used by the unconditionally 
!!!           stable and non-oscillatory BL numerical solver.  Also 
!!!           calculates screen level temperature and humidity as well 
!!!           as 10 m winds.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  6.4   10/01/07   New Deck         M. Diamantakis
!!!
!!!  Programming standard: UMDP4.
!!!
!!!
!!!  Documentation: 
!!!          http://www-nwp/~frmd/DR/Reports/new_BLsolver_guide.ps
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE SF_IMPL2 (                                             &

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
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,                         &
     & T_SOIL,QW_1,TL_1,U_1,V_1,RHOKM_U_1,RHOKM_V_1,GAMMA,              &
     & GAMMA1,GAMMA2,ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,               &

     & DTRDZ_CHARNEY_GRID_1,DU_1,DV_1,                                  &
     & FQW_TILE,EPOT_TILE,FQW_ICE,FTL_ICE,                              &
     & FRACA,RESFS,RESFT,RHOKH,RHOKH_TILE,RHOKH_SICE,                   &
     & RHOKPM,RHOKPM_POT,RHOKPM_SICE,Z1,                                &
     & Z0HSSI,Z0MSSI,Z0H_TILE,Z0M_TILE,CDR10M_U,CDR10M_V,               &
     & CHR1P5M,CHR1P5M_SICE,CT_CTQ_1,CTCTQ1,DQW_1,DTL_1,                &
     & DQW1_1,DTL1_1,DU_STAR1,DV_STAR1,CQ_CM_U_1,CQ_CM_V_1,             &
     & L_NEG_TSTAR,L_correct,FLANDG_U,FLANDG_V,ANTHROP_HEAT,            &
     & L_SICE_HEATFLUX,                                                 &

! IN STASH flags :-
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,                            &

! INOUT data :
     & TI,TI_GB,TSTAR,                                                  &
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,                       &
     & TSTAR_TILE,SNOW_TILE,                                            &
     & LE_TILE,RADNET_SICE,RADNET_TILE,                                 &
     & E_SEA,FQW_1,FTL_1,FTL_TILE,H_SEA,OLR,TAUX_1,TAUY_1,              &
     & TAUX_LAND,TAUX_LAND_star,TAUX_SSI,TAUX_SSI_star,TAUY_LAND,       &
     & TAUY_LAND_star,TAUY_SSI,TAUY_SSI_star,                           &

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
#if defined(ACCESS)
      USE auscom_cpl_data_mod,                                          &
     &      Only : auscom_salinity, ocn_sss, access_tfs
#endif

      IMPLICIT NONE

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
     &,NICE                        ! IN Number of sea ice catagories

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
     &,RAD_SICE(ROW_LENGTH,ROWS)                                        &
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
     & GAMMA1(row_length,rows)                                          &
                                   ! weights for new BL solver          
     &,GAMMA2(row_length,rows)                                          

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
     &,DU_1(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)                &
                                   ! IN Level 1 increment to u wind
!                                  !    field
     &,DV_1(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:N_ROWS+OFF_Y)              &
!                                  ! IN Level 1 increment to v wind
!                                  !    field
     &,ANTHROP_HEAT(NTILES)                                             &
                                   ! IN Additional heat source on tiles
                                   !    for anthropogenic urban heat
                                   !    source  (W/m2) 
     &,CTCTQ1(ROW_LENGTH,ROWS)                                          &
     &,DQW1_1(ROW_LENGTH,ROWS)                                          &
     &,DTL1_1(ROW_LENGTH,ROWS)                                          &
     &,DU_STAR1(1-off_x:row_length+off_x,1-off_y:rows+off_y)            &
     &,DV_STAR1(1-off_x:row_length+off_x,1-off_y:n_rows+off_y)          &
     &,TAUX_LAND_star(ROW_LENGTH,ROWS)                                  &
     &,TAUX_SSI_star(ROW_LENGTH,ROWS)                                   &
     &,TAUY_LAND_star(ROW_LENGTH,N_ROWS)                                &
     &,TAUY_SSI_star(ROW_LENGTH,N_ROWS)
!                                  ! IN Additional arrays needed by the 
!                                  !    uncond stable BL numerical solver


      LOGICAL                                                           &
     & LTIMER                                                           &
                                   ! IN Logical switch for TIMER diags
     &,L_NEG_TSTAR                                                      &
                                   ! IN Switch for -ve TSTAR error check
     &,L_FLUX_BC                                                        &
!                                  ! IN SCM logical for prescribed
                                   !    surface flux forcing
     &,L_correct                                                        &
                                   ! flag used by the new BL solver
     &,L_SICE_HEATFLUX             ! IN T: semi-implicit sea ice temp

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

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL IM_SF_PT2,SF_EVAP,SF_MELT,SCREEN_TQ,SICE_HTF
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

#include "c_0_dg_c.h"
#include "c_r_cp.h"
#include "c_lheat.h"
#include "csigma.h"
#include "caoptr.h"
#include "c_mdi.h"

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

      REAL ltfs

#if defined(ACCESS)
      ltfs = access_tfs
#else
      ltfs = tfs
#endif

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_IMPL2 ',3)
      ENDIF
      ERROR = 0

! DEPENDS ON: im_sf_pt2
      CALL IM_SF_PT2 (                                                  &
     &   OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,LAND_PTS                    &
     &,  LAND_INDEX,NTILES,TILE_INDEX,TILE_PTS                          &
     &,  FLANDG,TILE_FRAC,SNOW_TILE,ICE_FRACT                           &
     &,  GAMMA,GAMMA1,GAMMA2,ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE        &
     &,  RESFT,RHOKPM,RHOKPM_POT,RHOKPM_SICE                            &
     &,  RHOKM_U_1,RHOKM_V_1,RHOKH_TILE,RHOKH_SICE                      &
     &,  CT_CTQ_1,CTCTQ1,DQW_1,DTL_1,DQW1_1,DTL1_1                      &
     &,  CQ_CM_U_1,CQ_CM_V_1,DU_1,DV_1,DU_STAR1,DV_STAR1                &
     &,  FLANDG_U,FLANDG_V,FQW_1,FTL_1                                  &
     &,  TAUX_1,TAUX_LAND,TAUX_LAND_star,TAUX_SSI,TAUX_SSI_star,TAUY_1  &
     &,  TAUY_LAND,TAUY_LAND_star,TAUY_SSI,TAUY_SSI_star                & 
     &,  FQW_TILE,EPOT_TILE,FTL_TILE,FQW_ICE,FTL_ICE,E_SEA,H_SEA        &
     &,  L_CORRECT,L_FLUX_BC,LTIMER                                     &
     &  )
!
!-----------------------------------------------------------------------
!
! Calculate surface scalar fluxes, temperatures only at the 1st call
! of the subroutine (first stage of the new BL solver) using standard
! MOSES2 physics and equations. These are the final values for this
! timestep and there is no need to repeat the calculation.
!-----------------------------------------------------------------------
!
      IF ( .NOT. L_correct ) THEN
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

!-----------------------------------------------------------------------
! Diagnose the GBM surface temperature for points with sea-ice
!-----------------------------------------------------------------------

        TSTAR_SIC(:,:,:)= 0.0
        TSTAR_SSI0(:,:) = TSTAR_SSI(:,:)
        SURF_HT_FLUX_SICE(:,:)=0.0

        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          IF ( FLANDG(I,J) <  1.0 .and. ICE_FRACT(I,J) >  0. ) THEN
           SURF_HT_FLUX_SICE(I,J) = RADNET_SICE(I,J)-LS*FQW_ICE(I,J) -  &
     &                              FTL_ICE(I,J)
           TSTAR_SSI0(I,J) = (1. - ICE_FRACT(I,J))*TSTAR_SEA(I,J) +     &
     &                        ICE_FRACT(I,J)*TI(I,J,1) +                &
     &                        SURF_HT_FLUX_SICE(I,J) / ASHTF(I,J)
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

        DO N=1,NTILES
          DO L=1,LAND_PTS
            TSTAR_TILE(L,N) = T_SOIL(L,1)
            IF (SNOW_TILE(L,N) >  0.)                                   &
     &      TSTAR_TILE(L,N) =  MIN( T_SOIL(L,1), TM )
          ENDDO
        ENDDO

        DO N=1,NTILES
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            LAT_HT = LC
            IF (SNOW_TILE(L,N) >  0.) LAT_HT = LS
            TSTAR_TILE(L,N) = T_SOIL(L,1) + ( RADNET_TILE(L,N)          &
     &                          - LAT_HT*FQW_TILE(L,N) - FTL_TILE(L,N)  &
     &                          + (CANHC_TILE(L,N)/TIMESTEP) *          &
     &                            (TSTAR_TILE_OLD(L,N) - T_SOIL(L,1)) ) &
     &                                                 / ASHTF_TILE(L,N)
          ENDDO
        ENDDO

!-----------------------------------------------------------------------
!!     Surface evaporation components and updating of surface
!!     temperature (P245, routine SF_EVAP).
!-----------------------------------------------------------------------
! DEPENDS ON: sf_evap
        CALL SF_EVAP (                                                  &
     &   ROW_LENGTH,ROWS,LAND_PTS,NTILES,                               &
     &   LAND_INDEX,TILE_INDEX,TILE_PTS,SM_LEVELS,LTIMER,FLAND,         &
     &   ASHTF_TILE,CANOPY,DTRDZ_CHARNEY_GRID_1,FLAKE,FRACA,SNOW_TILE,  &
     &   RESFS,RESFT,RHOKH_TILE,TILE_FRAC,SMC,WT_EXT_TILE,TIMESTEP,     &
     &   FQW_1,FQW_TILE,FTL_1,FTL_TILE,TSTAR_TILE,                      &
     &   ECAN,ECAN_TILE,ELAKE_TILE,ESOIL,ESOIL_TILE,EI_TILE,EXT         &
     & )

!-----------------------------------------------------------------------
!!     Surface melting of sea-ice and snow on land tiles.
!-----------------------------------------------------------------------
! DEPENDS ON: sf_melt
        CALL SF_MELT (                                                  &
     &   ROW_LENGTH,ROWS,LAND_PTS,NTILES,NICE,LAND_INDEX,               &
     &   TILE_INDEX,TILE_PTS,LTIMER,SIMLT,SMLT,FLANDG,                  &
     &   ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,DTRDZ_CHARNEY_GRID_1,      &
     &   ICE_FRACT,RESFT,RHOKH_TILE,RHOKH_SICE,TILE_FRAC,TIMESTEP,      &
     &   GAMMA,                                                         &
     &   EI_TILE,FQW_1,FQW_ICE,FTL_1,FTL_ICE,FTL_TILE,                  &
     &   TSTAR_SEA,TSTAR_SIC,TSTAR_TILE,SNOW_TILE,                      &
     &   EI_LAND,EI_SICE,SICE_MELT,ICE_FRACT_NCAT,                      &
     &   SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,MELT_TILE,               &
     &   TSTAR_SSI0,SICE_MELT0 )

        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            QIM_1(I,J)=QW_1(I,J)+DQW1_1(I,J)-CTCTQ1(I,J)*FQW_1(I,J)
            TIM_1(I,J)=TL_1(I,J)+DTL1_1(I,J)-CTCTQ1(I,J)*FTL_1(I,J)/CP
            IF ( FLANDG(I,J) < 1.0 .and. ICE_FRACT(I,J) > 0. ) THEN
              IF(NICE == 1)THEN
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

!-----------------------------------------------------------------------
!!     Specific humidity and temperature at 1.5 metres.
!-----------------------------------------------------------------------
! DEPENDS ON: screen_tq
        CALL SCREEN_TQ (                                                &
     &   ROW_LENGTH,ROWS,LAND_PTS,NTILES,                               &
     &   LAND_INDEX,TILE_INDEX,TILE_PTS,FLANDG,                         &
     &   SQ1P5,ST1P5,CHR1P5M,CHR1P5M_SICE,PSTAR,QIM_1,RESFT,            &
     &   TILE_FRAC,TIM_1,TSTAR_SSI,TSTAR_TILE,                          &
     &   Z0HSSI,Z0H_TILE,Z0MSSI,Z0M_TILE,Z1,                            &
     &   Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE,                             &
     &   lq_mix_bl                                                      &
     & )

!-----------------------------------------------------------------------
!!     Gridbox-mean surface temperature and net surface heat fluxes
!-----------------------------------------------------------------------
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          SURF_HT_FLUX_LAND(I,J) = 0.
          SURF_HT_FLUX_SICE(I,J) = 0.
          DO N=1,NICE
            SURF_HT_FLUX_SICE_NCAT(I,J,N)=0.
          ENDDO
          TSTAR_RAD4(I,J) = 0.
          IF (FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0.) THEN
            TSTAR_SICE(I,J) = (TSTAR_SSI(I,J) -                         &
     &              (1.-ICE_FRACT(I,J))*TSTAR_SEA(I,J))/ICE_FRACT(I,J)
            TSTAR_RAD4(I,J) = (1.0-FLANDG(I,J))                         &
     &              *ICE_FRACT(I,J)*TSTAR_SICE(I,J)**4
            RADNET_SICE(I,J) = RAD_SICE(I,J) -                          &
     &                     ICE_FRACT(I,J)*SBCON*TSTAR_SICE(I,J)**4
            IF ( NICE == 1 ) THEN
              SURF_HT_FLUX_SICE(I,J) = RADNET_SICE(I,J)-LS*FQW_ICE(I,J)-&
     &                        FTL_ICE(I,J) - SICE_MELT0(I,J)
              SURF_HT_FLUX_SICE_NCAT(I,J,1) = SURF_HT_FLUX_SICE(I,J)
            ELSE
              DO N=1,NICE
                SURF_HT_FLUX_SICE_NCAT(I,J,N) = (ICE_FRACT_NCAT(I,J,N)/ &
     &                        ICE_FRACT(I,J)) * (RADNET_SICE(I,J) -     &
     &                        LS*FQW_ICE(I,J) -                         &
     &                        FTL_ICE(I,J)) - SICE_MELT(I,J,N)          &
     &                        + 4.0*SBCON*(TSTAR_SICE(I,J)**3.0) *      &
     &                        (ICE_FRACT_NCAT(I,J,N)*TSTAR_SICE(I,J)    &
     &                        - TSTAR_SIC(I,J,N))
                SURF_HT_FLUX_SICE(I,J) = SURF_HT_FLUX_SICE(I,J) +       &
     &                                   SURF_HT_FLUX_SICE_NCAT(I,J,N)
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

       If (LRAD_EMIS_LAND_GEN) Then

! ANTHROP_HEAT(N) takes value 0.0 except when N=6 (urban) and  
! L_ANTHROP_HEAT_SRC=.true. 
 
         DO N=1,NTILES
           DO K=1,TILE_PTS(N)
             L = TILE_INDEX(K,N)
             J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
             I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
             RADNET_TILE(L,N) = SW_TILE(L,N) + RAD_EMIS_LAND_GEN*       &
     &                        (LW_DOWN(I,J) - SBCON*TSTAR_TILE(L,N)**4)
             LE_TILE(L,N) = LC*ECAN_TILE(L,N) + LC*ESOIL_TILE(L,N) +    &
     &                      LC*ELAKE_TILE(L,N) + LS*EI_TILE(L,N)
             SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J)            &
     &                      + TILE_FRAC(L,N) *                          &
     &                      ( RADNET_TILE(L,N) + ANTHROP_HEAT(N) -      &
     &                        FTL_TILE(L,N) -                           &
     &                        LE_TILE(L,N) - LF*MELT_TILE(L,N) -        &
     &                       (CANHC_TILE(L,N)/TIMESTEP) *               &
     &                       (TSTAR_TILE(L,N) - TSTAR_TILE_OLD(L,N)) )
             TSTAR_LAND(I,J) = TSTAR_LAND(I,J)                          &
     &                + TILE_FRAC(L,N)*TSTAR_TILE(L,N)
             TSTAR_RAD4(I,J) = TSTAR_RAD4(I,J) + FLANDG(I,J)*           &
     &                     TILE_FRAC(L,N)*TSTAR_TILE(L,N)**4
           ENDDO
         ENDDO

       Else
        
         DO N=1,NTILES
           DO K=1,TILE_PTS(N)
             L = TILE_INDEX(K,N)
             J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
             I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
             RADNET_TILE(L,N) = SW_TILE(L,N) +                          &
     &                         LW_DOWN(I,J) - SBCON*TSTAR_TILE(L,N)**4
             LE_TILE(L,N) = LC*ECAN_TILE(L,N) + LC*ESOIL_TILE(L,N) +    &
     &                      LC*ELAKE_TILE(L,N) + LS*EI_TILE(L,N)
             SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J)            &
     &                      + TILE_FRAC(L,N) *                          &
     &                      ( RADNET_TILE(L,N) + ANTHROP_HEAT(N) -      &
     &                        FTL_TILE(L,N) -                           &
     &                        LE_TILE(L,N) - LF*MELT_TILE(L,N) -        &
     &                       (CANHC_TILE(L,N)/TIMESTEP) *               &
     &                       (TSTAR_TILE(L,N) - TSTAR_TILE_OLD(L,N)) )
             TSTAR_LAND(I,J) = TSTAR_LAND(I,J)                          &
     &                + TILE_FRAC(L,N)*TSTAR_TILE(L,N)
             TSTAR_RAD4(I,J) = TSTAR_RAD4(I,J) + FLANDG(I,J)*           &
     &                     TILE_FRAC(L,N)*TSTAR_TILE(L,N)**4
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
        CALL SICE_HTF(                                                  &
     &   ROW_LENGTH,ROWS,FLANDG,SIMLT,NICE,                             &
     &   DI_NCAT,ICE_FRACT,ICE_FRACT_NCAT,SURF_HT_FLUX_SICE_NCAT,       &
     &   TSTAR_SICE,TIMESTEP,TI,SICE_MLT_HTF,SEA_ICE_HTF,               &
     &   L_SICE_HEATFLUX,LTIMER )

!-----------------------------------------------------------------------
!  Calculate surface latent heat flux.
!-----------------------------------------------------------------------

        IF (SLH) THEN
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              LATENT_HEAT(I,J) = LC*FQW_1(I,J)                          &
     &     + LF*(FLANDG(I,J)*EI_LAND(I,J)+(1.-FLANDG(I,J))*EI_SICE(I,J))
            ENDDO
          ENDDO
        ENDIF

!-----------------------------------------------------------------------
!    Set RHOKH, the coefficients required for tracer mixing.
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

        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            TSTAR(I,J)=FLANDG(I,J)*TSTAR_LAND(I,J)                      &
     &        +(1.-FLANDG(I,J))*TSTAR_SSI(I,J)
            EI(I,J)=FLANDG(I,J)*EI_LAND(I,J)                            &
     &        +(1.-FLANDG(I,J))*EI_SICE(I,J)

            SURF_HT_FLUX(I,J)=FLANDG(I,J)*SURF_HT_FLUX_LAND(I,J)        &
     &        +(1.-FLANDG(I,J))*SURF_HT_FLUX_SICE(I,J)

            IF ( NICE  > 1 )THEN
              IF (FLANDG(I,J) == 1.0) THEN
                TI_GB(I,J) = RMDI
              ELSE IF (ICE_FRACT(I,J) <= 0.0) THEN
#if defined(ACCESS)
                if (ocn_sss) then
                    ltfs = ZeroDegC - 0.054 * auscom_salinity(I,J)
                end if
#endif
                TI_GB(I,J) = lTFS
              ELSE
                TI_GB(I,J) = 0.0
                DO N=1,NICE
                  TI_GB(I,J)=TI_GB(I,J)+ICE_FRACT_NCAT(I,J,N)           &
     &                      /ICE_FRACT(I,J)*TI(I,J,N)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO

!-----------------------------------------------------------------------
! Rescale FTL_1 as it should be used to update the botom row of the 
! discrete equation handled by the new BL solver at the next (2nd)
! stage of the scheme. 
!-----------------------------------------------------------------------
!fpp$ Select(CONCUR)
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            FTL_1(I,J) = FTL_1(I,J)/CP
          ENDDO
        ENDDO

      ELSE ! L_correct = true: 2nd stage of the scheme

!-----------------------------------------------------------------------
! Rescale to Watts/m^2 as this is the final call to the imp BL solver
! and FTL_1 will be used by stash diagnostics
!-----------------------------------------------------------------------
!fpp$ Select(CONCUR)
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            FTL_1(I,J) = CP*FTL_1(I,J)
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
!  U_V will be updated at 2nd stage of the scheme as the equations
!  providing the implicit surface stresses have been modified
!  consistently with the new scheme.
!-----------------------------------------------------------------------
! U component of 10m wind
        IF (SU10) THEN
          DO J=1,ROWS
           DO I=1,ROW_LENGTH
            U10M(I,J) = (U_1(I,J) + DU_STAR1(I,J) + (DU_1(I,J) -        &
     &                   CQ_CM_U_1(I,J)*TAUX_1(I,J)) -                  &
     &                   U_0(I,J))*CDR10M_U(I,J) + U_0(I,J)             
           ENDDO
          ENDDO
        ENDIF

! V component of 10m wind
        IF (SV10) THEN
          DO J=1,N_ROWS
            DO I=1,ROW_LENGTH
              V10M(I,J) = (V_1(I,J) + DV_STAR1(I,J) + (DV_1(I,J) -      &
     &                     CQ_CM_V_1(I,J)*TAUY_1(I,J)) -                &
     &                     V_0(I,J))*CDR10M_V(I,J) + V_0(I,J)
            ENDDO
          ENDDO
        ENDIF

! Correct surface stress diagnostics

        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            TAUX_LAND(I,J) = TAUX_LAND(I,J) + TAUX_LAND_star(I,J)
            TAUX_SSI(I,J)  = TAUX_SSI(I,J)  + TAUX_SSI_star(I,J)
          ENDDO
        ENDDO

        DO J=1,N_ROWS
          DO I=1,ROW_LENGTH
            TAUY_LAND(I,J) = TAUY_LAND(I,J) + TAUY_LAND_star(I,J)
            TAUY_SSI(I,J)  = TAUY_SSI(I,J)  + TAUY_SSI_star(I,J)
          ENDDO
        ENDDO

      END IF ! IF .NOT. L_correct

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_IMPL2 ',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_IMPL2
#endif
